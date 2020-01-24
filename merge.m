clear all;
close all;
files= dir('csi*.dat');
lf=length(files);
data = cell(1,lf);
for i=1:lf
    filename=files(i).name;
    data{i}=importdata(filename);
end
tot_data=zeros(size(data{1}));
for i=1:lf
    tot_data=tot_data+data{i};
end
tot_data=tot_data/lf;
dlmwrite('csi.dat',tot_data, 'delimiter','\t');

exit;