% Simulation of spot commodity prices. Use for Oultook Forecast.
% Mean reversion decay factor pass into as input parameter
% Run this separately for each fuel, changing the names of the output file.

function [Returned_70thDist Returned_Shock] = PriceForecastOutlookCorr (output_name, name, ...
    spot_price, volatility, MeanRevRate, MeanRevDecayFactor, Seed, randoms)

global dates;
global num_runs;
global days_in_month;
global rep_monthly;
global rep_annual;
global num_days

num_months=length(spot_price);  %number of months in the spotprice set
num_years = ceil(num_months/12);%number of years in the spotprice set 
MeanRevScalingFactor=1;

% Call function that reads input spot_price array
% and expand the monthly data to the daily set
forecast = expandForecast(spot_price);

% copy expanded dataset 
exp_spot_price = forecast; 
long_term_means = forecast;

start_price = exp_spot_price(1);

% create place holder for random generated numbers
%randoms=zeros(num_runs,num_days);
t=zeros(num_days+1,1);
MeanReversion=zeros(num_days,1);

%start clock to monitor execution time
fix(clock);
tic;

distribution=zeros(num_runs,num_days+1);
distribution(:,1)=start_price;
t(1)=1;

% % % % % % FirstDay_of_Month = xlsread('Reset_Monthly_MeanValues.xlsx','FirstDay_of_Month','A2:A266');
count=1;

for day=1:num_days
    t(day+1)=day+1;

% % % % % %     if (day == FirstDay_of_Month(count,1))
% % % % % %         %we build the shock in parts... start with the reversion to the mean
% % % % % %         shock=(log(long_term_means(day))-log(long_term_means(day)))*MeanRevRate; %start with reversion%
% % % % % %     else
        shock=(log(long_term_means(day))-log(distribution(:,day)))*MeanRevRate; %start with reversion%
% % % % % %     end
    
    %smoothing decay line
    MeanReversion(day) = Seed/(day+Seed)*(1-MeanRevDecayFactor)+MeanRevDecayFactor;
    
    %linear growth of reversion strength     
    shock=shock*MeanReversion(day);
    
    %...finally, add the volatility-scaled random shock
    shock=shock+randoms(:,day)*volatility;
    tempshock(:,day) = shock;
    
% % % % % %     if (day == FirstDay_of_Month(count,1))
% % % % % %         distribution(:,day) = spot_price(count);
% % % % % %         distribution(:,day+1)=(distribution(:,day).*exp(shock));%+drift_adder;
% % % % % %         count=count+1;
% % % % % %     else
        distribution(:,day+1)=(distribution(:,day).*exp(shock));%+drift_adder;        
% % % % % %     end
end
    Returned_Shock = tempshock(:,:);

% Get percentiles of the distribution
prctiles=prctile(distribution, [1 5 25 50 75 95 99]);
Returned_70thDist = distribution;

% Combine bands with forecast prices
CL_bands=[exp_spot_price prctiles'];
   
disp([num2str(sum(CL_bands(:,1)>CL_bands(:,7))/num_days*100) '%  above 95%']);
disp([num2str(sum(CL_bands(:,1)<CL_bands(:,3))/num_days*100) '%  below 5%']);

% Monthly distribution
tic;
monthly_distribution=zeros(num_runs,(num_days+1)/days_in_month);
for trail=1:num_runs  
  for counter=1:(num_days+1)/days_in_month
      monthly_distribution(trail,counter)= sum(distribution(trail,(counter-1)*days_in_month+1:counter*days_in_month))/days_in_month;
  end
end
  
% Get percentiles of distribution
mnthly_prctiles=prctile(monthly_distribution, [1 5 25 50 75 95 99]);
month_CL_bands=[spot_price mnthly_prctiles'];

% % % % % % disp([num2str(sum(month_CL_bands(:,1)>month_CL_bands(:,7))/(num_days+1)/days_in_month*100) '%  above 95%']);
% % % % % % disp([num2str(sum(month_CL_bands(:,1)<month_CL_bands(:,3))/(num_days+1)/days_in_month*100) '%  below 5%']);
% % % % % % col_header={'Forecast ', 'Low 1% ', 'Low 5% ', 'Low 25% ', 'Expected ', 'High 75% ', 'High 95% ', 'High 99% '}; %Row cell array (for column labels)

if (rep_monthly == 1)
% % % % % %     % Combine headers, titles and save simulated distribution to Excel
% % % % % %     %
% % % % % %     month_CL_bands=[ {' '} col_header; dates num2cell(month_CL_bands)]; %Join cell arrays
% % % % % %     
% % % % % %     disp('Writing monthly bands to Excel');
% % % % % %     file_name=sprintf('Monthly_%s.xls',output_name );  
% % % % % %     warning off MATLAB:xlswrite:AddSheet;
% % % % % %     xlswrite(file_name, month_CL_bands, name);
% % % % % %     pause(6);
end

% Annual distribution
if (rep_annual == 1)
    [year month] = datevec(dates);
    ind = find(month==12,1);

    annual_dist=zeros(num_runs,ceil(num_months/12));

    % Get an average forecast for each month
    if (ind ~=12) % if first year is not a full 12 months
        avg_forecast(1)=sum(spot_price(1:ind))/ind;
        for j=1:length(unique(year))-1 
            avg_forecast(j+1)=sum(spot_price((j-1)*12+1+ind:min(j*12+ind,length(month))))/(min(j*12+ind,length(month))-((j-1)*12+ind));
        end    
    else
         for j=1:length(unique(year))         
              avg_forecast(j)=sum(spot_price((j-1)*12+1:min(j*12,length(month))))/(min(j*12,length(month))-((j-1)*12));
         end 
    end

    % Get an average prices for each month
    for i=1:num_runs
         if (ind ~=12) % if first year is not a full 12 months
             annual_dist(i,1)= sum(monthly_distribution(i,1:ind))/ind;
             
                for j=1:length(unique(year))-1                    
                    annual_dist(i,j+1)=sum(monthly_distribution(i,(j-1)*12+1+ind:min(j*12+ind,length(month))))/(min(j*12+ind,length(month))-((j-1)*12+ind));
                end    
         else
                for j=1:length(unique(year))        
                   annual_dist(i,j)= sum(monthly_distribution(i,(j-1)*12+1:min(j*12,length(month))))/(min(j*12,length(month))-((j-1)*12));
                end   
         end 
    end   

    % Get percentiles of distribution
    annual_prctiles=prctile(annual_dist, [1 5 25 50 75 95 99]);
    annual_CL_bands=[avg_forecast' annual_prctiles'];
%     Combine headers, titles and save simulated distribution to Excel
%     annual_CL_bands=[ {' '} col_header; num2cell(unique(year)) num2cell(annual_CL_bands)]; %Join cell arrays
%     disp('Writing annual bands to Excel');
%     file_name=sprintf('Annual_%s.xls',output_name );  
%     xlswrite(file_name, annual_CL_bands, name);
%     pause(6);
end

% This function reads input spot_price array and expand the monthly
% data to the daily by repetition of spot_price for each day of the months.
function x = expandForecast(spot_price)
global days_in_month;

num_months=length(spot_price);
forecast=zeros(days_in_month*num_months,1);
month=ones(days_in_month,1);

for counter=1:num_months
   forecast((counter-1)*days_in_month+1:counter*days_in_month)=month*spot_price(counter);
end

x = forecast;