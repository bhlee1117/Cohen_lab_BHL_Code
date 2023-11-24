function PlotSpeciesComparison( S , cutoff,listnames)
%PlotSpeciesComparison
%Plots the species data in the first two column of data versus each other.
%Adds names to the plots based on the cutoff of required diffrence between
%the values.
%If listnames is not 0 list names on plot

figure,plot(S.data(:,1),S.data(:,2),'.','Color',0.75*[1 1 1])
%The array following color is the amount of [red green blue] in the color
%The numbers can be between 0 and 1.  All 0 is black.  All 1 is white.
hold on
SignificantORFs=abs(S.data(:,1)-S.data(:,2))>cutoff;
%In this case the array data is assumed to be in log scale already
plot(S.data(SignificantORFs,1),S.data(SignificantORFs,2),'.','Color',0.25*[1 1 1])
if listnames~=0
    text(S.data(SignificantORFs,1),S.data(SignificantORFs,2),S.ORF(SignificantORFs))
    xlabel(S.header(1));
    ylabel(S.header(2));
end

end

