% re-ordering legends for an open plot
%% legends 
clear all
clc

% get the handle of the current plots
 h = get(gca,'Children');
 
 % Then select handles of lines you want to add to legend.
 % Say h(1) is handle of a line in 1st group,
 %    h(10)                       2nd group,
 %    h(15)                       3rd
 %    h(30)                       4th 
 v1 = [ h(1) h(6) h(5) h(4) h(7) h(2)  h(3)  ]';
 legend(v1);
