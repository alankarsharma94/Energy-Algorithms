%function SpotSimCorr

global dates;
global num_runs;            %number of simulations to run
global days_in_month;       %trading days, not calendar days
global rep_monthly;
global rep_annual;
global num_days

uiwait(msgbox('This program propagates spot prices (single or multiple commodities) and returns monthly and annual percentile bands','modal'));
%Input file
uiwait(msgbox('Select input file','modal'));
[filename,pathname] = uigetfile('*.xlsx','Input file');
input_file_name = fullfile(pathname,filename);

uiwait(msgbox('Select input parameters','modal'));
prompt={'Number of runs','Days in month','Monthly report (1:Yes, 0:No)','Annual report (1:Yes, 0:No)','Output file name'};

name='Parameters for spot price propagation';
numlines=1;
defaultanswer={'1000','30','1','0','2018Q3_Update'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);

num_runs = str2double(answer(1));
days_in_month = str2double(answer(2));
rep_monthly = str2double(answer(3));
rep_annual = str2double(answer(4));
output_name = char(answer(5));

disp('Starting Spot simulation');

volatility = xlsread(input_file_name, 'Volatility');
MeanRevRate = xlsread(input_file_name, 'Mean reversion');
Corr = xlsread(input_file_name, 'Correlation');
[forecast timeseries] = xlsread(input_file_name, 'Forecast');
dates = timeseries(2:end,1);
name = timeseries(1,2:end);

% Choleski decomposition
[LTDecom err]= chol(Corr);

while err ~=0
    Corr = Rebonato(Corr);
    [LTDecom err]= chol(Corr);    
end

num_days = size(forecast,1)*days_in_month-1;
randMx = zeros(length(name), num_days*num_runs);

for i=1:num_runs
    %populate randoms matrix with uncorrelated numbers 
    clear randoms;
    randoms=normrnd(0,1,[num_days length(name)]);
% % % %     %**********************************************************************
% % % %     randoms=randn([num_days length(name)]);
% % % %     %Upper and Lower limits for the random numbers
% % % %     random_lowerlimit = -1;
% % % %     random_upperlimit = 1;
% % % %     randoms(randoms<random_lowerlimit) = random_lowerlimit + (0.05-0.01)* randn([1 1])+0.01;
% % % %     randoms(randoms>random_upperlimit) = random_upperlimit + (0.05-0.01)* randn([1 1])+0.01;
% % % %     %**********************************************************************
    
    %transform uncorrelated random number to correlated using
    %Cholesky decomposition 
    rdCor = randoms*LTDecom;
    randMx(:,num_days*(i-1)+1:num_days*(i-1)+num_days) = rdCor';
    %clear rdCor;
end    
 
 
% % % for i = 1:length(volatility)
% % %     commodity_name = char(name(i));
% % %     randomsC = reshape(randMx(i,:), num_days, num_runs); %correlated random numbers per commodity
% % %     randomsC = randomsC';
% % % 
% % %     disp(['Working on ' commodity_name]); 
% % %     Returned_70thDist = PriceForecastOutlookCorr(output_name, commodity_name, forecast(:,i), volatility(i), MeanRevRate(1,i), MeanRevRate(2,i), MeanRevRate(3,i),randomsC);
% % % end
    
for i = 1:3  
    randomsC = reshape(randMx(i,:), num_days, num_runs); %correlated random numbers per commodity
    randomsC = randomsC';

    switch i
        case 1
            %HH Ref Case
            commodity_name = 'HHRef_70thPctage';
            disp(['Working on ' commodity_name]); 
            [Returned_70thDist_HHRef MainRet_Shock] = PriceForecastOutlookCorr('HHRef_70thPctage', ...
                'HHRef', forecast(:,i), volatility(1), MeanRevRate(1,1), MeanRevRate(2,1), MeanRevRate(3,1),randomsC);
            %CaseRefShock = MainRet_Shock;
        case 2
            %HH Low Case
            commodity_name = 'HHLow_15thPctage';
            disp(['Working on ' commodity_name]); 
            [Returned_Low15thDist_HHLow MainRet_Shock] = PriceForecastOutlookCorr('HHLow_15thPctage', ...
                'HHLow15', forecast(:,i), volatility(1), MeanRevRate(1,1), MeanRevRate(2,1), MeanRevRate(3,1),randomsC);
            %CaseLowShock = MainRet_Shock;
        case 3
            %HH High Case
            commodity_name = 'HHHigh_15thPctage';
            disp(['Working on ' commodity_name]); 
            [Returned_Upper15thDist_HHHigh MainRet_Shock] = PriceForecastOutlookCorr('HHHigh_15thPctage', ...
                'HHUpper15', forecast(:,i), volatility(1), MeanRevRate(1,1), MeanRevRate(2,1), MeanRevRate(3,1),randomsC);
            %CaseHighShock = MainRet_Shock;
    end
end
disp('Spot simulation has finished');

for GasPoints = 1:1  
    switch GasPoints
        case 1
            titlename = 'HH';
            OPFilename = 'Output_HH_Monthly_Dist_Transposed.xlsx';
            Returned_70thDist = Returned_70thDist_HHRef;
            Returned_Low15thDist = Returned_Low15thDist_HHLow;
            Returned_Upper15thDist = Returned_Upper15thDist_HHHigh;
    end
            
    Refcase_Permutations = randperm(1000);  %1000
    Lowcase_Permutations = randperm(150); %150
    Highcase_Permutations = randperm(150); %150

    ColIndex_Refcase = Refcase_Permutations(1,1:700); %1:700
    ColIndex_Lowcase = Lowcase_Permutations(1,1:150); %1:150
    ColIndex_Highcase = Highcase_Permutations(1,1:150); %1:150

    FinalMatrix=zeros(num_runs,(num_days+1));
    for i=1:size(ColIndex_Refcase,2)
        FinalMatrix(i,:) = Returned_70thDist(ColIndex_Refcase(1,i),:);
    end

    a=1;
    for j = (i+1):(i+size(ColIndex_Lowcase,2))
        FinalMatrix(j,:) = Returned_Low15thDist(ColIndex_Lowcase(1,a),:);
        a=a+1;
    end

    b=1;
    for k = (j+1):(j+size(ColIndex_Highcase,2))
        FinalMatrix(k,:) = Returned_Upper15thDist(ColIndex_Highcase(1,b),:);
        b=b+1;
    end

    %Convert Daily values to Monthly
    monthly_distribution=zeros(num_runs,(num_days+1)/days_in_month);
    for trail=1:700  %1:700  
      for counter=1:(num_days+1)/days_in_month
          monthly_distribution(trail,counter)= sum(FinalMatrix(trail,(counter-1)*days_in_month+1:counter*days_in_month))/days_in_month;
      end
    end

    for trail=701:850  %701:850
      for counter=1:(num_days+1)/days_in_month
          monthly_distribution(trail,counter)= sum(FinalMatrix(trail,(counter-1)*days_in_month+1:counter*days_in_month))/days_in_month;
      end
    end

    for trail=851:1000  %851:1000
      for counter=1:(num_days+1)/days_in_month
          monthly_distribution(trail,counter)= sum(FinalMatrix(trail,(counter-1)*days_in_month+1:counter*days_in_month))/days_in_month;
      end
    end

    %Randomize the data from the 3 cases (i.e. 700 iterations from RefCase,
    %150 iterations from LowCase and 150 iterations from HighCase)
    %They should not occur in a sequence in the final matrix.
    Rand_Sequence = randperm(1000);  %randperm(1000)
    for permutationcount = 1:size(Rand_Sequence,2)
        Final_Monthly_Dist(permutationcount,:) = monthly_distribution(Rand_Sequence(1,permutationcount),:);
    end
    xls(OPFilename,Final_Monthly_Dist,'Monthly_Prices','B2');
    
    
    %Impose floor on monthly prices not to go below $2.75 for HH
    if (GasPoints==1)
        Final_Monthly_Dist(Final_Monthly_Dist>15.00)=(15 - ((3-1)*randn(1,1)+0.1));
    else
        Final_Monthly_Dist(Final_Monthly_Dis<2)=(2 + ((0.3-0.1)*randn(1,1)+1));
    end

    %PDF Plot
    %Monthly Prices Distribution
    rn1_Series = Final_Monthly_Dist(:);
    murn1=mean(rn1_Series(:,1));
    sigmarn1=std(rn1_Series(:,1));
    Pdf_rn1=normpdf(rn1_Series(:,1),murn1,sigmarn1);

    if (GasPoints==1)
        figure;plot(rn1_Series(:,1),Pdf_rn1,'r*');    
    else
        hold;plot(rn1_Series(:,1),Pdf_rn1,'b*');
    end
    legend('HH Prices');

    %Convert Monthly to Yearly
    %[1000 X 276]
    for irows = 1:num_runs
        for icols = 1:(size(Final_Monthly_Dist,2)/12)
            Yearly_Price(irows,icols)= mean(Final_Monthly_Dist(irows,(((icols-1)*12+1):((icols-1)*12+12))));
        end
    end
    %monthly_distribution_trans = Final_Monthly_Dist';
    %xls(OPFilename,Yearly_Price,'Yearly_Prices','B2');

 end

 disp('Final Matrix Populated');