clc;
clear;
s=4 ; %Number of neurons
 
%Train Datas
train_data=-10:0.5:10; %Train Input Data
train_data=train_data.';
y_real_train=sin(train_data)./train_data; %Output of Train Input Data from System
y_real_train(21)=1; %To prevent the error at train input data equals to 0.


[X,y_model_value]=SISO_ANN_Train(s,train_data,y_real_train); %Training Process


%Test Datas
test_data=-8.5:0.75:8.75; %Test Input Data
test_data=test_data.';
y_real_test=sin(test_data)./test_data; %Output of Test Input Data from System



[y_m]=SISO_ANN_Test(X,s,test_data,y_real_test); %Test Process
