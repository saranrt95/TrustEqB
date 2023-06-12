clear all
clc
%load the excel file with the number of hits
data=readtable("NumberOfHits_fi34_fi05.xlsx");
%% Plot Histograms

for c=1:size(data,2)
    hit=data{:,c};
    figure
    bar(hit);
    xlabel("Rules");
    xticks([1:size(data,1)]);
    ylabel("Number of Hits");
    axis tight
    if c<=5
        title("Training domain");
    else
        title("Operational domain");
    end
    saveas(gcf, strcat("bar_plots_fi01_nuovefi34_basenuove/hitrate_ID",num2str(c),".png"));
end

%% Baseline range
% indexes of the first Histogram
ix34=[1,1,1,1,2,2,2,3,3,4];
% indexes of the second Histogram
iy34=[2,3,4,5,3,4,5,4,5,5];

for i=1:length(ix34)
        hitx=data{:,ix34(i)};
        hity=data{:,iy34(i)};
        %mi34(i)=mutualInformation(hitx,hity);
        mi34(i)=mutual_info(hitx,hity);
end


disp("Training baseline")
min(mi34)
max(mi34)


%% Operational range with approach 1

for c=1:5
    ix=[1,1,1,1,1];
    ix=c.*ix;
    iy=[6,7,8,9,10];
    for i=1:length(ix)
            hitx=data{:,ix(i)};
            hity=data{:,iy(i)};
            mi(i)=mutual_info(hitx,hity);
            MImatrix(i,c)=mi(i);
    end
    
    
end

disp("Operational baseline")
min(MImatrix(:))
max(MImatrix(:))

%% Operational range with approach 2

ixOP=[6,6,6,6,7,7,7,8,8,9];
iyOP=[7,8,9,10,8,9,10,9,10,10];

for i=1:length(ixOP)
        hitx=data{:,ixOP(i)};
        hity=data{:,iyOP(i)};
        miOP(i)=mutual_info(hitx,hity);
end
disp("Operational baseline")
min(miOP)
max(miOP)


%% Plot mutual information with approach 1
idx = 5;
matrixMI = MImatrix(:);
% da qua possiamo visualizzare se siamo OOD
figure
p1=plot(MImatrix(:),'LineWidth',1);
hold on
p2=yline(min(mi34),'color',[0 0.5 0],'LineWidth',1);
hold on
p3=yline(max(mi34),'r','LineWidth',1);
hold on
yline(min(MImatrix(:)),'k--','LineWidth',1);
hold on
yline(max(MImatrix(:)),'m--','LineWidth',1);
xticks([1:25]);
xlabel("Histograms couple",'FontSize',12,'FontWeight','bold');
ylabel("Mutual Information (\mu I)",'FontSize',12,'FontWeight','bold');
l= legend('$\mathbf{\mu I}_{op}$', '$min(\mathbf{\mu I}_{tr})$', '$max(\mathbf{\mu I}_{tr})$', '$min(\mathbf{\mu I}_{op})$', '$max(\mathbf{\mu I}_{op})$','FontSize',14,'FontWeight','bold');
l.Interpreter='latex';

%% Simple statistics (mean, variance, skewness and kurtosis)

dataMat=cell2mat(table2cell(data));
mean_hits=mean(dataMat);
var_hits=var(dataMat);
skew_hits=skewness(dataMat);
kurt_hits=kurtosis(dataMat);
hitfile=data.Properties.VariableNames;
stats=table(hitfile',mean_hits',var_hits',skew_hits',kurt_hits','VariableNames',["HitFileID","Mean","Variance","Skewness","Kurtosis"]);


%% Plot statistics

figure 
subplot(2,2,1)
plot(mean_hits,'LineWidth',1);
hold on
plot(1:5,mean_hits(1:5),'b*');
hold on
yline(mean(mean_hits(1:5)),'Color','k','LineStyle','--','LineWidth',1,'HandleVisibility','off');
hold on
plot(6:10,mean_hits(6:10),'r*');
hold on 
yline(mean(mean_hits(6:10)),'Color','k','LineStyle','-.','LineWidth',1,'HandleVisibility','off');
xlabel("Test iteration");
xlim([1,10])
ylabel("Mean");
l=legend("$m(\mathbf{w}^j$)","$\overline{\phi}_{34}$ domain","$\overline{\phi}_{01}$ domain","FontSize",14);
l.Interpreter='latex';

subplot(2,2,2)
plot(var_hits,'LineWidth',1);
hold on
plot(1:5,var_hits(1:5),'b*');
hold on
yline(mean(var_hits(1:5)),'Color','k','LineStyle','--','LineWidth',1,'HandleVisibility','off');
hold on
plot(6:10,var_hits(6:10),'r*');
hold on 
yline(mean(var_hits(6:10)),'Color','k','LineStyle','-.','LineWidth',1,'HandleVisibility','off');
xlabel("Test iteration");
xlim([1,10])
ylabel("Variance")
l=legend("$\sigma^2(\mathbf{w}^j$)","$\overline{\phi}_{34}$ domain","$\overline{\phi}_{01}$ domain","FontSize",14);
l.Interpreter='latex';

subplot(2,2,3)
plot(skew_hits,'LineWidth',1);
hold on
plot(1:5,skew_hits(1:5),'b*');
hold on
yline(mean(skew_hits(1:5)),'Color','k','LineStyle','--','LineWidth',1,'HandleVisibility','off');
hold on
plot(6:10,skew_hits(6:10),'r*');
hold on 
yline(mean(skew_hits(6:10)),'Color','k','LineStyle','-.','LineWidth',1,'HandleVisibility','off');
xlabel("Test iteration");
xlim([1,10])
ylabel("Skewness")
l=legend("$s(\mathbf{w}^j$)","$\overline{\phi}_{34}$ domain","$\overline{\phi}_{01}$ domain","FontSize",14);
l.Interpreter='latex';

subplot(2,2,4);
plot(kurt_hits,'LineWidth',1);
hold on
plot(1:5,kurt_hits(1:5),'b*');
hold on
yline(mean(kurt_hits(1:5)),'Color','k','LineStyle','--','LineWidth',1,'HandleVisibility','off');
hold on
plot(6:10,kurt_hits(6:10),'r*');
hold on 
yline(mean(kurt_hits(6:10)),'Color','k','LineStyle','-.','LineWidth',1,'HandleVisibility','off');
xlabel("Test iteration");
xlim([1,10])
ylabel("Kurtosis")
l=legend("$\kappa(\mathbf{w}^j$)","$\overline{\phi}_{34}$ domain","$\overline{\phi}_{01}$ domain","FontSize",14);
l.Interpreter='latex';

