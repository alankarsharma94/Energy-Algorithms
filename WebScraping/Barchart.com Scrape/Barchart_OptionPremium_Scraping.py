# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 08:20:27 2020

@author: Alankar Sharma
"""
import requests
import codecs
import json
import pandas as pd
from lxml import html
from urllib.parse import unquote

#This function sends the HTTP request in order to grab the first expiration option URL from Barchart
def grab_first_expiration(site):
    #Send the GET request
    qq = requests.get(site, headers=headers, stream=False, timeout=(8, 15))

    # Check if the cookie in the headers is empty, this implies that the first request is to grab the cookie
    if headers['Cookie'] == '':
        cookie = ''
        # Grab all the cookies from the response
        for x in qq.cookies.get_dict():
            cookie += x + '=' + qq.cookies.get_dict()[x] + ';'
        # Set the XSRF-TOKEN token, this is required for future requests to work
        headers['x-xsrf-token'] = unquote(qq.cookies.get_dict()['XSRF-TOKEN'])
        # Set the cookie from above
        headers['Cookie'] = cookie

    #Set the response's text into an HTML string, we are going to parse it below
    tree = html.fromstring(qq.text)

    #Grab all the first expiration date as URL
    return tree.xpath("//div[@class='error-page']/ul/li/a/@href")[0]

#This function sends the HTTP request in order to grab the data from Barchart
def response(site, ticker, expire, readjson):
    #Send the GET request
    qq = requests.get(site, headers=headers, stream=False, timeout=(8, 15))

    if readjson:
        #If the cookie has been set, then we are just grabbing the data. The data is a JSON format
        json_check(qq.text, ticker, expire)
    else:
        # Set the response's text into an HTML string, we are going to parse it below
        tree = html.fromstring(qq.text)

        # Grab all the Option expiration dates as URLs, so we can loop them
        for value in tree.xpath("//select[@id='bc-options-toolbar__dropdown-month']//option/@value"):
            l.append(value)

        # Grab all the Option expiration dates as text, so we can loop them
        for value in tree.xpath("//select[@id='bc-options-toolbar__dropdown-month']//option/text()"):
            l2.append(value)

#This function creates the CSV file
def write_file(string_to_write): 
    with codecs.open(r'O:\PTIFairfax\Risk Management\Options Data\GlobOptions.csv', mode='w', encoding='utf-8') as f: #Update where you want to save
        f.write(string_to_write)

#This function parses the JSON from the response function
def json_check(json_value, ticker, expire):
    #This loads the JSON
    x = json.loads(json_value)

    global csv_to_write

    #Loop the JSON response and add it the CSV string
    for options in x['data']:
        csv_to_write += ''.join([i for i in options['strike'].replace(',','') if i.isdigit() or i == '.']) + ','
        csv_to_write += options['highPrice'].replace(',','') + ','
        csv_to_write += options['lowPrice'].replace(',','') + ','
        csv_to_write += ''.join([i for i in options['lastPrice'].replace(',','') if i.isdigit() or i == '.']) + ','
        csv_to_write += options['priceChange'].replace(',','') + ','
        csv_to_write += options['bidPrice'].replace(',','') + ','
        csv_to_write += options['askPrice'].replace(',','') + ','
        csv_to_write += options['volume'].replace(',','') + ','
        csv_to_write += options['openInterest'].replace(',','') + ','
        csv_to_write += options['premium'].replace(',','') + ','
        csv_to_write += options['tradeTime'].replace(',','') + ','
        csv_to_write += ''.join([i for i in options['strike'].replace(',', '') if not (i.isdigit() or i == '.')]) + ','
        csv_to_write += ''.join([i for i in options['lastPrice'].replace(',', '') if not (i.isdigit() or i == '.')]) + ','
        csv_to_write += expire.replace(' ',' 1 ') + '\n'

#These are headers to be sent with HTTP requests. For now, this is just the User-Agent header, to make the request appear as if it's coming from a browser.
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.88 Safari/537.36',
    'Cookie': '',
}

l = list()
l2 = list()

#This creates the CSV headers
csv_to_write = 'strike,highPrice,lowPrice,lastPrice,priceChange,' \
   'bidPrice,askPrice,volume,openInterest,premium,tradeTime,StrikeIdentifier,lastPriceIndentifier,contractDate' + '\n'

url = 'https://www.barchart.com' + grab_first_expiration('https://www.barchart.com/futures/quotes/NGJ20/options/JLNJ20')

#Send the first request. This is just to grab the cookie and the Option expiration dates.
response(url, '', '', False)

#Loop each Option expiration dates.
for list_value, list_value2 in zip(l, l2):
    #Add these headers. The Referer is not required. The Accept is required in order to receive the JSON response.
    headers['Accept'] = 'application/json'
    headers['Referer'] = 'https://www.barchart.com' + list_value

    #Parse the ticker and expiration date
    ticker = list_value[-6:]
    expire = list_value2

#####Send the request to Barchart, and then we'll parse the response##########
#######I clicked Ctrl + Shift + I (do this before the page loads, or reload the page after this step) > clicked the network tab > I clicked Ctrl + F to pull up the find feature >
############ I entered one of the option values (1.4500 for instance) > only one search came up, so I clicked that link
    response('https://www.barchart.com/proxies/core-api/v1/quotes/get?symbol={0}&list=futures.options&fields=strike%2ChighPrice%2ClowPrice%2ClastPrice%2CpriceChange%2CbidPrice%2CaskPrice%2Cvolume%2CopenInterest%2Cpremium%2CtradeTime%2ClongSymbol%2CsymbolCode%2CsymbolType%2ChasOptions&meta=field.shortName%2Cfield.description%2Cfield.shortName%2Cfield.type&hasOptions=true&raw=1'.format(ticker), ticker, expire, True)

write_file(csv_to_write)

######################### Getting Final Data ##################################
#GlobOptions.csv created manually for the first round of data pull to initialize dataframe and column names

Data = pd.read_csv(r"GlobOptions.csv") #### This calls the path that you saved GlobOptions make sure to update it

#Getting the first 5 rows and listing all columns names
Data.head(5)
for col in Data.columns:
    print(col)

#Renaming the data and dropping some variables
InProcessData = pd.DataFrame(Data)
InProcessData.drop(["priceChange", "bidPrice", "askPrice", "volume", "openInterest", "premium"], axis = 1, inplace = True)

#Transforming the contractDate to MM-DD-YYYY
ContractDate = pd.to_datetime(InProcessData['contractDate'])
print(ContractDate)

#Adding ContractDate to the data and dropping the old contractDate
InProcessData['ContractDate'] = ContractDate
print(InProcessData)

InProcessData.drop(['contractDate'], axis = 1, inplace = True)
for col in InProcessData.columns:
    print(col)
InProcessData.head(5)

#Getting the sum and average of high and low price (use avg when last price not available)
InProcessData['Sum_High_Low'] = InProcessData.highPrice.add(InProcessData.lowPrice, fill_value=0)
InProcessData['Average_high_Low'] = (InProcessData['Sum_High_Low'])/2

#Replacing missing from lastPrice to the value of corresponding row of Average_high_Low 
InProcessData.lastPrice = InProcessData.lastPrice.fillna(InProcessData.Average_high_Low)

#Dropping highPrice, lowPrice, and Average_high_Low_Price
InProcessData.drop(["highPrice", "lowPrice", "Sum_High_Low", "Average_high_Low"], axis = 1, inplace = True)

#Saving as csv file and naming it FinalData
FinalData = pd.DataFrame(InProcessData)
FinalData.to_csv(r'FinalData_Barchart.csv', index=False) #######Update the path