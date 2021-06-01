# Flexible Resource Adequacy (Subhourly Timescale) Estimation Model
# V2
# Author A. Sharma
# Date 9/8/2020

## Methodology ##

# Read 1 Hour Load, Read renewable capacities, Hourly Profiles
# Sub-Hourly standard deviations for each resource type is calculated externally through pivot tables
# Sub-Hourly Deviations applied to hourly profiles to create normally distributed variations at 1 min level
# Cubic interpolation applied between hour ending timestamps to also create sub-hourly level profile
# Resource Adequacy requirement is the percentile of difference of sub-hourly variation and smoothed sub-hourly profile

import sys
import os
import os.path

from datetime import datetime
import numpy
import random
import numpy as np
from numpy import array
import random
from scipy import interpolate
from xlrd import xldate_as_tuple
import pandas as pd
from openpyxl import load_workbook
import logging
import h5py
import time

#Fit a cubic spline function to interpolate minute to minute variation into hourly load
def Trans1hourto1minNoVar(y):
    nrows = len(y)
    x = np.linspace(1,nrows,nrows)
    f = interpolate.interp1d(x,y,'cubic')
    xnew = np.linspace(1,(nrows),(nrows)*60)
    ynew = f(xnew)
    return ynew


iFile = 0 # don't change this - used to read in sections of large files
Iter = 5 #change total Iterations to evaluate here. Moer iterations, longer runtimes


def runFlex(Iter, iFile):
    InputFile = 'RE Inputs.xlsx' #change input file name correspondingly
    iSection = iFile * Iter # Legacy - dont change it

    totEntries = 19 ;  # change number of plants plus load
    totCase = 6;  # change total case number
    totHour = 8784; # change total hour if leap year 8784 if non-leap 8760
    totMinute = totHour*6;
    MaxPercentile = 95; # change percentile to evaluate
    MinPercentile = 100 - MaxPercentile;

    """
    Also change 1-min peak load around line 159 to properly account for max load and load factors
 
    """


    #----------------------------------------
    startTime = datetime.now()
    logger = logging.getLogger()
    handler = logging.FileHandler('RA_Log' + str(iFile)+'.txt') #save the error logging file to the same path
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    logger.info('\n ------------------ Reporting tool started on %s ----------------------' %str(startTime))
    logger.info('\n')

    try:

        idxx = 0
        idRE = 0
        idHr = 0
        CaseId = 0

        df_Shape = pd.read_excel(InputFile,sheet_name = 'HourlyShape')
        df_Shape.dropna(how='all')
        df_Vara = pd.read_excel(InputFile,sheet_name = 'HourlyVariation')
        df_1min = pd.read_excel(InputFile,sheet_name = 'OneMinTime')
        df_Cap = pd.read_excel(InputFile,sheet_name = 'RE Capacity')

        df_hrindex = pd.read_excel(InputFile,sheet_name = 'HourIndex',usecols="A:C")
        df_intervalindex = pd.read_excel(InputFile,sheet_name = '10minIndex',usecols="A:D")
        #Create backups in pickle format
        df_hrindex.to_pickle('df_hrindex.pkl')
        df_intervalindex.to_pickle('df_intervalindex.pkl')
        df_Shape.to_pickle('df_Shape.pkl')
        df_Vara.to_pickle('df_Vara.pkl')
        df_1min.to_pickle('df_1min.pkl')

        RE_Capacity = df_Cap[df_Cap.columns[1:totCase+1]].values # column 0 - 5 represents case 1 - 6
        RECapacity=list(zip(*RE_Capacity))
        RECapacity=np.asarray(RECapacity)

        np.save('RECapacity.npy',RECapacity)

        RE_Hourly = []
        for i in range(0,totEntries):
            RE_Hourly.append(df_Shape[df_Shape.columns[i+3]].values)

        RE_1min = []
        for i in range(0,totEntries):
            RE_1min.append(Trans1hourto1minNoVar(RE_Hourly[i])) #cubic interpolate hourly data to 1-min data

        end = str(datetime.now())
        logger.info('----------------- 1min data creation ended on %s ----------------------' %end)

        np.save('REHourly.npy', RE_Hourly)
        np.save('RE1min.npy',RE_1min)


        """
        Load and 15 Resources including DG ;  16 in total
        """

        errmin = 0
        errmax = 0
        err_RE = 0

        """
        #REIter[] 1st layer - 1,000 iterations data - RE_EachIter
        #         2nd layer - 16 Entries _RE_EachSite
        #         3rd layer - 8784*60 minutes data
        """

        REIter = [] # initialize 1st layer

        Varaidx = [] #varaibility index
        for j in range(0,totHour): # 8760 if not leap year
            iMonth = df_Shape.iloc[j,1]
            iHour = df_Shape.iloc[j,2]

            idx = df_Vara.loc[(df_Vara['MonthIndex'] == iMonth) & (df_Vara['HourIndex'] == iHour)].index
            Varaidx.append(idx)

        for i in range(0, Iter): #1st layer 5 iteration
            idxx = i
            RE_EachIter = []
            for m in range(0,totEntries): #2nd layer 16 entries;
                idRE = m
                RE_EachSite = []
                for j in range(0,totHour): # 3rd layer 8784 hours
                    idHr = j
                    if (RE_1min[m][j*60] < 0.0001) and (RE_1min[m][59+j*60] < 0.0001): #get rid of extreme low number in raw inputs
                        for n in range(0,60):
                            RE_EachSite.append(0)
                    else: #otherwise, get normal distribution by applying mean and std. dev.

                        idx = Varaidx[j]
                        sigma_RE = df_Vara.iloc[idx,m+2].values[0]

                        for n in range(0,60): #4th layer 60 minutes
                            RE1min = RE_1min[m][n+j*60]

                            if m == 0:
                                errmax = min(3289/RECapacity[0][0],(RE1min + 3*sigma_RE)) # change 1-min peak load; 2976 is the 1-min peak load provided by client
                                #errmax = min(1,(RE1min + 3*sigma_RE)) # if no 1-min peak load defined, use this line and comment out the line above
                            else:
                                errmax = min(1,(RE1min + 3*sigma_RE))
                            errmin = max(0,(RE1min - 3*sigma_RE))

                            err_RE = random.normalvariate(RE1min, sigma_RE)

                            err_RE = max(err_RE,errmin)
                            err_RE = min(err_RE,errmax)
                            err_RE = max(err_RE, 0)
                            RE_EachSite.append(err_RE)

                RE_EachIter.append(RE_EachSite)

                del RE_EachSite

            """
            change below to plug in the section if needed
            """
            with h5py.File('REIter'+str(i+iSection)+'.h5', 'w') as hf:
                hf.create_dataset('REIter'+str(i+iSection),  data=array(RE_EachIter))

            del RE_EachIter


        del Varaidx
        del RE1min

        time.sleep(120)

        """
        Only exporting below
        """
##        df_hrindex = pd.read_pickle('df_hrindex.pkl') #If need to recall from backup
##        df_intervalindex= pd.read_pickle('df_intervalindex.pkl') #If need to recall from backup
##
##        df_1min= pd.read_pickle('df_1min.pkl') #If need to recall from backup

        hpArr = []
        for iter in range(0,Iter*(iFile+1)):
            with h5py.File('REIter'+str(iter)+'.h5', 'r') as hf:
                Arr = []
                Arr = hf['REIter'+str(iter)][:]
            hpArr.append(Arr)
            del Arr
        PerCase = []
        minNLArr = []
        for CaseId in range(0,totCase): # 6 is the total case number (6 capacity scenarios to evaluate)
            minNL = []
            RECap = RECapacity[CaseId]
            for iter in range (0,Iter*(iFile+1)):

                NetLoad = hpArr[iter][0]*(RECap[0])

                for m in range(1,totEntries):
                    NetLoad = NetLoad - hpArr[iter][m]*RECap[m]

                df = pd.DataFrame(data=list(NetLoad),columns=['NetLoad'])
                dfNL = df_1min.join(df)


                df_HourlyNL = dfNL.pivot_table(index='HourofYear',values = 'NetLoad',aggfunc = [np.mean])
                df_10minNL = dfNL.pivot_table(index=['HourofYear','10min Interval'],values = 'NetLoad', aggfunc = [np.mean])

                if iter == 0:
                    df_hrNlIter = df_HourlyNL.copy()
                    df_10minNLIter = df_10minNL.copy()
                else:
                    df_hrNlIter = df_hrNlIter.join(df_HourlyNL,rsuffix=str(iter))
                    df_10minNLIter = df_10minNLIter.join(df_10minNL,rsuffix=str(iter))

                tenminNL = np.asarray(df_10minNL.values)

                HourlyNL = np.asarray(df_HourlyNL.values)


                HourlyNLDelta = HourlyNL[3:]-HourlyNL[:-3] #calculate three-hour net load delta
                HourlyNLDelta =np.append(HourlyNLDelta,0)# to get 8784 record
                HourlyNLDelta =np.append(HourlyNLDelta,0)# to get 8784 record
                HourlyNLDelta =np.append(HourlyNLDelta,0)# to get 8784 record

                tenminDelta = tenminNL[1:]-tenminNL[:-1] #52703 record
                tenminDelta = np.append(tenminDelta,0) # to get 52704 record

                df_hrdelta = pd.DataFrame(data=list(HourlyNLDelta),columns=['HourlyDelta'])
                df_10mindelta = pd.DataFrame(data=list(tenminDelta),columns=['10minDelta'])

                if iter == 0:
                    df_hrdeltaIter = df_hrdelta.copy()
                    df_10mindeltaIter = df_10mindelta.copy()
                else:
                    df_hrdeltaIter = df_hrdeltaIter.join(df_hrdelta,rsuffix=str(iter))
                    df_10mindeltaIter = df_10mindeltaIter.join(df_10mindelta,rsuffix=str(iter))


                del NetLoad
                del df_HourlyNL, df_10minNL, df_10mindelta,df_hrdelta,HourlyNLDelta, HourlyNL, tenminDelta, tenminNL, m

            """
            Below is to output all iterations of 10 min net load data
            """

            df_10minNLIter.to_pickle('10minNLIter_'+str(CaseId)+'.pkl')
            df_hrNlIter.to_pickle('hrNLIter_'+str(CaseId)+'.pkl')
            df_10minNLIter.to_csv('output_10minNL_'+str(CaseId)+'.csv') ## can comment out to save time


            df_hrNlIter.reset_index(inplace=True)


            df_hrdeltaIter = df_hrindex.join(df_hrdeltaIter)
            df_10mindeltaIter = df_intervalindex.join(df_10mindeltaIter)
            df_hrNlIter = df_hrindex.join(df_hrNlIter) # for max and min



            df_hrNlIter = df_hrNlIter[:totHour]
            df_hrdeltaIter = df_hrdeltaIter[:totHour]
            df_10mindeltaIter = df_10mindeltaIter[:totMinute]

            """
            Below is to output all iterations of hourly net load data
            """

            df_hrNlIter.to_csv('output_hrNL_'+str(CaseId)+'.csv') ## can comment out to save time

            """
            Below is to output all iterations of hourly net load DELTA data
            """

            df_hrdeltaIter.to_csv('output_hrdelta_'+str(CaseId)+'.csv')
            df_10mindeltaIter.to_csv('output_10mindelta_'+str(CaseId)+'.csv')

            df_hrNlIter1 = df_hrNlIter.query('Month>4 & Month < 10') # select month May to Sep change month if needed
            df_hrNlIter2 = df_hrNlIter.query('(Month<5 & (Month > 1)) | (Month ==10) | (Month == 11)') #select Month Feb-Apri OR Oct. OR Nov. change month if needed

            nar = df_hrNlIter1[df_hrNlIter1.columns[4:]].to_numpy(copy=False)
            nar2 = df_hrNlIter2[df_hrNlIter2.columns[4:]].to_numpy(copy=False)
            nar3 = df_hrNlIter2[df_hrNlIter2.columns[4:]].to_numpy(copy=False)


            maxnumber=numpy.percentile(nar,MaxPercentile,interpolation = 'nearest')
            minnumber=numpy.percentile(nar2,MinPercentile,interpolation = 'nearest')
            minNL = numpy.sum(nar3[np.where(nar3<0)])/1000
            irow, jrow = np.where(df_hrNlIter.values == maxnumber)
            irow2, jrow2 = np.where(df_hrNlIter.values == minnumber)

            RAresults = []
            RAresults.append(df_hrNlIter.loc[irow[0],'HourlyIndex'])
            RAresults.append(df_hrNlIter.loc[irow[0],'Month'].item(0))
            RAresults.append(df_hrNlIter.loc[irow[0],'Hour'])
            RAresults.append(maxnumber)


            RAresults.append(df_hrNlIter.loc[irow2[0],'HourlyIndex'])
            RAresults.append(df_hrNlIter.loc[irow2[0],'Month'].item(0))
            RAresults.append(df_hrNlIter.loc[irow2[0],'Hour'])
            RAresults.append(minnumber)

            del nar, nar2, maxnumber, minnumber,irow,irow2, jrow, jrow2
            del df_hrNlIter1, df_hrNlIter2, df_hrNlIter


            for imonth in range(0,12):

                df_hrdeltaIter1 = df_hrdeltaIter[df_hrdeltaIter['Month'] == imonth+1].copy()
                nar = df_hrdeltaIter1[df_hrdeltaIter1.columns[3:]].to_numpy(copy=False)
                #Calc monthly percentile levels from subhourly data
                maxnumber=numpy.percentile(nar[np.where(nar>0)],MaxPercentile,interpolation = 'nearest')
                minnumber=numpy.percentile(nar[np.where(nar<0)],MinPercentile,interpolation = 'nearest')
                irow, jrow = np.where(df_hrdeltaIter.values == maxnumber)
                irow2, jrow2 = np.where(df_hrdeltaIter.values == minnumber)


                RAresults.append(df_hrdeltaIter.loc[irow[0],'HourlyIndex'])
                RAresults.append(df_hrdeltaIter.loc[irow[0],'Month'].item(0))
                RAresults.append(df_hrdeltaIter.loc[irow[0],'Hour'])
                RAresults.append(maxnumber)


                RAresults.append(df_hrdeltaIter.loc[irow2[0],'HourlyIndex'])
                RAresults.append(df_hrdeltaIter.loc[irow2[0],'Month'].item(0))
                RAresults.append(df_hrdeltaIter.loc[irow2[0],'Hour'])
                RAresults.append(minnumber)

                del nar, maxnumber, minnumber,irow,irow2, jrow, jrow2
                del df_hrdeltaIter1

                df_10mindeltaIter1 = df_10mindeltaIter[df_10mindeltaIter['Month'] == imonth+1].copy()
                nar = df_10mindeltaIter1[df_10mindeltaIter1.columns[4:]].to_numpy(copy=False)
                maxnumber=numpy.percentile(nar[np.where(nar>0)],MaxPercentile,interpolation = 'nearest')
                minnumber=numpy.percentile(nar[np.where(nar<0)],MinPercentile,interpolation = 'nearest')

                irow, jrow = np.where(df_10mindeltaIter.values == maxnumber)
                irow2, jrow2 = np.where(df_10mindeltaIter.values == minnumber)

                RAresults.append(df_10mindeltaIter.loc[irow[0],'10minInterval'])
                RAresults.append(df_10mindeltaIter.loc[irow[0],'Month'].item(0))
                RAresults.append(df_10mindeltaIter.loc[irow[0],'Hour'])
                RAresults.append(maxnumber)
                RAresults.append(df_10mindeltaIter.loc[irow2[0],'10minInterval'])
                RAresults.append(df_10mindeltaIter.loc[irow2[0],'Month'].item(0))
                RAresults.append(df_10mindeltaIter.loc[irow2[0],'Hour'])
                RAresults.append(minnumber)


                del nar
                del df_10mindeltaIter1, irow,irow2, jrow, jrow2, maxnumber, minnumber

            PerCase.append(RAresults)
            minNLArr.append(minNL)
            del RAresults
            del minNL
            del df_10mindeltaIter, df_hrdeltaIter

        del hpArr, RECap

        dfFinal = pd.DataFrame(PerCase,columns=['Timestamp','Month','Hour','MaxSeasonalNetLoad',
                                                'Timestamp','Month','Hour','MinSeasonalNetLoad',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn1',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn1',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn2',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn2',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn3',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn3',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn4',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn4',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn5',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn5',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn6',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn6',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn7',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn7',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn8',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn8',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn9',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn9',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn10',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn10',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn11',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn11',
                                                'Timestamp','Month','Hour','Max3hourNetLoad_mn12',
                                                'Timestamp','Month','Hour','Min3hourNetLoad_mn12',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn1',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn1',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn2',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn2',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn3',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn3',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn4',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn4',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn5',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn5',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn6',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn6',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn7',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn7',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn8',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn8',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn9',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn9',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn10',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn10',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn11',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn11',
                                                'Timestamp','Month','Hour','Max10minNetLoad_mn12',
                                                'Timestamp','Month','Hour','Min10minNetLoad_mn12'])

        dfFinal.T.to_excel('Stat_p99.xlsx')

        df_minNL = pd.DataFrame(minNLArr,columns=['minimumNetLoadSum'])
        df_minNL.T.to_excel('MinNetLoadSum_p99.xlsx') # sum of all negative net load

        end = str(datetime.now())
        logger.info('----------------- Reporting tool ended on %s ----------------------' %end)
        logger.removeHandler(handler)
        print('success')

    except Exception as e:
        print(f'Error at case {CaseId} iter {idxx} RE {idRE} Timestamp {idHr} ')
        print(str(e))
        print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
        if (e != 0):
            logger.info('Error has been found when genearting RE data. Please check following error message.')
            logger.info(str(e))
            end = str(datetime.now())
            logger.info('----------------- Reporting tool ended on %s ----------------------' %end)
            logger.removeHandler(handler)


runFlex(Iter,iFile)
