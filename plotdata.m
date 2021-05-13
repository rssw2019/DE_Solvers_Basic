fileID = fopen('fort.200','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
A = fscanf(fileID,formatSpec,sizeA);
plot(A(1,:),A(2,:))

fclose(fileID);
