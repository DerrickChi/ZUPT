function assignment2(file)
% file = 'BlueRadios_6243_1375383815730.txt';
data = ReadData(file);
i = DetermineInitialWinodw(data);
acc = preprocess(data, i);
data = [data acc];
d = trajectory_estimation(acc, data);
norm(d(end,:)-d(1,:))
end

% read data from file
function data = ReadData(file)    
as = 2048;       % sensitivity of 16g accelerometer
gs = 16.4;       % sensitivity of gyroscope
qs = 1073741824; % sensitivity of quaternion
fid = fopen(file,'r');
raw_data = [];
id = 0;
a = fgets(fid);
while(ischar(a))
    id = id + 1;
    if id == 1
        a = fgets(fid);
        continue
    else
        raw_data = [raw_data; str2num(a)];
    end
    a = fgets(fid);
end
fclose(fid);
if ~isempty(raw_data)
    raw_data(:,2) = -raw_data(:,2);
    raw_data(:,3) = -raw_data(:,3);
    acc = raw_data(:,2:4)/as;
    gyro = raw_data(:,5:7)/gs;
    q = raw_data(:,8:11)/qs;
    data = [acc gyro q];
else
    data = [];
end
end

% determine the window for initialization
function i = DetermineInitialWinodw(data)
wlen = 20;
int_threshold = 0.002; 
wlen = 2*wlen;
for i = 1:wlen+1:length(data(:,1))
    acc = data(i:i+wlen,:);
    acc_var = var(acc(:,1).^2 + acc(:,2).^2 + acc(:,3).^2);
    if (acc_var >= int_threshold)
        break;
    end
end
i = i - 1;
if i == 0
    i = 1;
else
    i = i - wlen/2;
end
acc = data(:,1:3);
gyro = data(:,4:6);
q = data(:,7:10);
figure;
subplot(3,1,1);
plot(acc);
hold on;
plot([i i], get(gca,'YLim'), '--k');
hold off;
title('Sensor Measurment');
ylabel('Accelerometer [g]');
legend('X-axis', 'Y-axis', 'Z-axis', 'Window for Initialization', 'location','NorthWest');
subplot(3,1,2);
plot(gyro);
hold on;
plot([i i], get(gca,'YLim'), '--k');
hold off;
ylabel('Gyroscope [deg/s]');
subplot(3,1,3);
plot(q);
hold on;
plot([i i], get(gca,'YLim'), '--k');
hold off;
ylabel('Orientation [Quaternion]');
xlabel('Sample Index');
end

% map acceleration to global coordinate and perform gravity subtraction
function acc = preprocess(data, wlen)
acc = data(1:wlen,1:3);
q = data(1:wlen,7:10);
a_ref = mean(acc,1);
q_ref = mean(q,1);
acc = data(:,1:3);
q = data(:,7:10);
q = quatdivide(q, q_ref);                             
acc = gravity_subtraction(acc, a_ref, q);
figure;
subplot(2,1,1);
plot(data(:,1:3));
title('Acceleration');
ylabel('Acceleration in Sensor Cooridnate');
subplot(2,1,2);
plot(acc);
ylabel('Acceleration After Gravity Subtraction in Global Coordinate');
xlabel('Sample Index');
end

% subtrac gravity
function global_acc = gravity_subtraction(acc, acc_mean, q)
global_acc = zeros(size(acc));
for i = 1:length(acc)
    a_temp = [0 acc(i,:)];
    a_temp = quatmultiply(quatmultiply(q(i,:),a_temp),quatconj(q(i,:)));
    global_acc(i,:) = a_temp(2:4) - acc_mean;
end
end

% estimate foot motion
function d = trajectory_estimation(acc, data)
dt = 1/200;
a = acc*9.8;
v = cumsum(a)*dt;
v1 = zero_velocity_update(data, v);
d = cumsum(v1)*dt;
figure;
subplot(3,1,1);
plot(a);
title('Position Estimation');
ylabel('Acceleration [m^2/s]');
subplot(3,1,2);
plot(v1);
ylabel('Velocity [m/s]');
subplot(3,1,3);
plot(d);
ylabel('Displacement [m]');
xlabel('Sample Index');
end

% code for zero velocity update
function v1 = zero_velocity_update(data, v)
% data(:,1:3) accelerometer measurements
% data(:,4:6) gyroscope measurements
% data(:,7:10) quaternion representing orientations
% data(:,11:13) acceleration after gravity correction

% figure;
% subplot(4,1,1);
% plot(data(:,1:3));
% subplot(4,1,2);
% plot(data(:,4:6));
% subplot(4,1,3);
% plot(data(:,7:10));
% subplot(4,1,4);
% plot(data(:,11:13));
v1 = v;
end