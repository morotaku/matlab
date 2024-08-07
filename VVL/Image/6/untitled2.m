imp1=importdata('total.csv');
Imp1=imp1./max(max(imp1));

imagesc(Imp1)