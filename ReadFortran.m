% 从Fortran写好的文件中读取数据

%xa = 0; xb = 2*pi; ya = 0; yb = 2*pi; % Tang Vortex
%xa = -10; xb = 10; ya = -10; yb = 10; % 平移
xa = 0;xb = 1;ya = 0;yb = 1; % 云
%xa = 0;xb = 2*pi;ya = 0;yb = 2*pi; % sin
%xa = 0;xb = 10;ya = 0;yb = 10; % 等熵涡

% Q1 = str2num(fileread('Q1.txt'));
% Q1 = Q1(:,1);
% Q2 = str2num(fileread('Q2.txt'));
% Q2 = Q2(:,1);
% Q3 = str2num(fileread('Q3.txt'));
% Q3 = Q3(:,1);
% Q4 = str2num(fileread('Q4.txt'));
% Q4 = Q4(:,1);
% T = fileread('T.txt');
% T = str2num(T);

Q1 = readmatrix('Q1');
Q1 = Q1(:,1);
Q2 = readmatrix('Q2');
Q2 = Q2(:,1);
Q3 = readmatrix('Q3');
Q3 = Q3(:,1);
Q4 = readmatrix('Q4');
Q4 = Q4(:,1);
T = str2num(fileread('T.txt'));

N = round(sqrt(length(Q1)/length(T)));
% 网格
h = (xb - xa)/N;
h1 = h/2;
X = zeros(N,1);
Y = zeros(N,1);
for i = 1:N + 1
    X(i) = xa + (i - 1)*h;
    Y(i) = ya + (i - 1)*h;
end

% 整格点
Xc = (X(1:end - 1) + X(2:end))/2;
Yc = (Y(1:end - 1) + Y(2:end))/2;

Q1 = reshape(Q1,N,N,length(T));
Q2 = reshape(Q2,N,N,length(T));
Q3 = reshape(Q3,N,N,length(T));
Q4 = reshape(Q4,N,N,length(T));

gamma = 1.4;

QP = (gamma - 1)*(Q4 - 0.5*(Q2.^2 + Q3.^2)./Q1);

FlashFortran