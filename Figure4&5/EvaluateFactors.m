function EvaluateFactors( data )

% Pull the Z-data ouf of the data cell
data=data.spec;
%run svd
[~,S,~] = svd(data);

%Get the diagonal elements of the S matrix
Sdiag = diag(S,0);

%calculate the relative values of the first ten diagonal elements of the S
%matrices
sdiag_rel = ones(1,10);
for i = 1:10
    sdiag_rel(i) = Sdiag(i) / Sdiag(i+1);
end

%calculate the percentage of the data accounted for by each factor
[Sdiag_rows, ~] = size(Sdiag);
Sdiag_percent = ones(Sdiag_rows,1);
for i = 1:Sdiag_rows
    Sdiag_percent(i) = Sdiag(i)/sum(Sdiag);
end

%Generate a plot of the relative values of the first ten diagonal elements
%of the S matrices
xaxis = 1:10;
clf
subplot(2,1,1)
plot(xaxis,sdiag_rel);
xlabel('Factor');
ylabel('Relative S Value');
xlim([1 10])

%Generate a plot of the percent of data accounted for by each factor
subplot(2,1,2)
plot(xaxis,Sdiag_percent(1:10),xaxis,cumsum(Sdiag_percent(1:10)));
xlabel('Factor');
ylabel('Fraction of Data Accounted For');
legend('Fraction per Factor','Cumulative Sum')
ylim([0 1])
xlim([1 10])
drawnow
end

