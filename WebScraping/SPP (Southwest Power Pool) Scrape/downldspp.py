# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:32:23 2021

@author: z00435su
"""

import pandas
import time
import os
#wget is the key
import wget 
import ssl

#name of subfolder to save files
dirnm = 'SPPDAT'


def spp_dl(year,leap,dnm):
    
    #datetime stuff
    months=list(range(1,13))
    hours= ['0100'] #For day-ahead-fuel-on-margin dataset, daily level data. Otherwise create list with hour range
    
    #create the loop to make urls and download csv/zip files
    yr=year
    for mo in months:
        if mo == 1 or mo==3 or mo==5 or mo==7 or mo==8 or mo==10 or mo==12:
            dayct = 31
        if mo == 4 or mo ==6 or mo==9 or mo==11:
            dayct = 30
        if mo ==2:
            if leap==0:
                dayct=28
            else:
                dayct=29
        days = list(range(1,dayct+1))
        for da in days:
            for hr in hours:
                mo2 = str(mo).zfill(2)
                da2 = str(da).zfill(2)
                hr2 = str(hr).zfill(4)
                #create url with date components flowing from loop earlier
                url = "https://marketplace.spp.org/file-browser-api/download/day-ahead-fuel-on-margin?path=%2F{yr}%2F{mo2}%2FDA-FUEL-ON-MARGIN-{yr}{mo2}{da2}{hr2}.csv".format(yr=yr,mo2=mo2,da2=da2,hr2=hr2)

                fnm = str(yr)+str(mo2)+str(da2)+str(hr2) + '.csv'
                fpath = os.path.join(os.path.dirname(__file__),r"{}/{}".format(dnm,fnm))


                #download the data file
                founddata=False #True if connection established, false if error occured
                filedir = os.path.dirname(fpath)
                if not os.path.exists(filedir):
                    os.makedirs(filedir) #Create file directory to download into
                xlimit = 1
                while founddata==False and xlimit<10:
                    try:
                        print(filedir)
                        ssl._create_default_https_context = ssl._create_unverified_context
                        wget.download(url,out=filedir)
                        founddata=True
                        break
                    except:
                        print("Connection refused") #MOst likely error with URL string if error caught here
                        print("Wait for a short while")
                        xlimit=xlimit+1
                        print(url)
                        time.sleep(3)
                        continue

#call to spp_dl
#argument 1 = year
#argument 2 = leap year (1-yes; 0-no)
#argument 3 = the name of the subdirectory to save downloads specified at the top
spp_dl(2020,1,dirnm) #Call 2020 with leap year flag
spp_dl(2021,0,dirnm) #Call 2021 without leap year flag

    
