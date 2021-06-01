# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 02:47:37 2021

@author: Alankar Sharma
"""
from selenium import webdriver #Main package for web scraping
from lxml import html #For querying HTML webpage address 
import time #For pausing between code sections
from selenium.webdriver.support.select import Select
from webdriver_manager.chrome import ChromeDriverManager
from webdriver_manager.utils import ChromeType
import pandas as pd
from tqdm import tqdm #For progress bar

# URL to query - Main CME frontend for Henry Hub Option data (ProductID 1352 = European options)
url="https://www.cmegroup.com/trading/energy/natural-gas/\
natural-gas_quotes_globex_options.html#optionProductId=1352&strikeRange=ALL"

#Create function to initialize webpage
def initWebdriver(mode):
    options = webdriver.ChromeOptions()
    if(mode):
        options.add_argument('headless') #Headless mode to not open Chrome window
    else:
    	print("Running in Normal Mode") #Normal mode opens up a Chrome window with current page being pulled
    #Configure google chrome
    driver = webdriver.Chrome(ChromeDriverManager(chrome_type=ChromeType.GOOGLE).install(),options=options) 
    return driver

#Manual column creation for organizing queried dataset
columns=['CALLS UPDATED','CALLS VOLUME','CALLS HIGH','CALL SLOW',\
'CALLS PRIOR SETTLE','CALLS CHANGE','CALLS LAST',\
 'PUTS LAST','PUTS CHANGE','PUTS PRIOR SETTLE',\
 'PUTS LOW','PUTS HIGH','PUTS VOLUME','PUTS UPDATED','STRIKE PRICE','Date']

#Create function to query for data table from the open webpage
def getthetable(index,date):
	final=[]
	ele=Select(driver.find_element_by_id('cmeOptionExpiration')) #Query based on Option Expiration Dates (prompt month to furthest available date)
	ele.select_by_index(index) #select_by_index to query from dropdown list
	tree=html.fromstring(driver.page_source) #Retrieve webpage data for currently open webpage
    
    #Loop to append each option expiration data table to last queried data table based on xpath under observation
	for idx,row in enumerate(tree.xpath('//*[@id="optionQuotesProductTable1"]/tbody/tr')):
		rowdata=[]
		for i in row.xpath('./td'):
		    rowdata.append(i.text_content())
		rowdata.append(''.join(row.xpath('./th/text()')))
		rowdata.append(date)
		final.append(rowdata)
	return final

#COnfigure drivers
driver=initWebdriver(False) # Mode = False for tracking current open page
driver.get(url)
time.sleep(5)

tree=html.fromstring(driver.page_source)

#Call functions and create final_df with all expiration date's option data
finaldf=pd.DataFrame()
for index,date in enumerate(tqdm(driver.find_element_by_id('cmeOptionExpiration').text.split('\n'))):
	try:
		data=getthetable(index,date)
		df=pd.DataFrame(data,columns=columns)

	except Exception as e:
		print("Error occured",e)
	time.sleep(1)
	finaldf=pd.concat([df,finaldf])

finaldf.to_csv("Output.csv",index=False)
