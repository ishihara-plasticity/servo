clear
format compact
close all

%q=[x,xd,θ,θd,β,y,yd,Φ,Φd,γ,z,zd]'
% パラメータ
r = 0.6; % =(2/3)(振子の長さ)
g = 9.81; % 重力加速度 [m/s^2]

A = [0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 g 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 g/r 0 -g/r 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 -g 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 g/r 0 -g/r 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0;];

B = [0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    1 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 1 0;
    0 0 0;
    0 0 1;];

C = [1 0 0 0 0 0 0 0 0 0 0 0 0];
Ct = [1 0 0 0 0 0 0 0 0 0 0 0 0];

At = [A ,zeros(12, 1);
      -C,   0       ];
Bt = [B; zeros(1, 3)];

%重みの決定
%        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
%        [x,xd, θ,θd, β, y,yd, Φ,Φd,  γ,  z, zd, 追従偏差w]
Q = diag([1, 1, 1, 100, 1, 1, 1, 1, 100, 1, 1, 1, 6000]);
R=diag([1,1,1]);

%lqr法
Kt = lqr(At, Bt, Q, R);
K1 = Kt(:, 1:12);
K2 = Kt(:, 13);

% シミュレーションの設定
Tf = 120;  % シミュレーション時間 [s]
Ts = 0.05; % サンプリング時間 [s]
T = 0:Ts:Tf;  % シミュレーションの時間軸
N = length(T);  % シミュレーションの時間ステップ数(=201)

% 目標値
r0 = 1.0;

% 初期状態などの設定
q0 = [0; 0; 0.1;  0; 0; 0; 0; -0.15; 0; 0; 0; 0; 0]; % 初期状態として、振子を傾けておく

% 使用する変数等を0で宣言
q = zeros(13, N);
dq = zeros(13, N);
q(:, 1) = q0;
u = zeros(3, N); % 制御入力の初期化
q_ref = zeros(13, N);
u1 = zeros(3,1);
u2 = zeros(3,1);
a = 0;
da = 0;
y = zeros(13, N);

for k = 1:N-1

    % 状態フィードバックによる制御入力の計算
    u1(:,k) = -K1*(q(1:12, k) - q_ref(1:12, k));
    integralDiff = Ts*q(13, k);
    u2(:,k) = -K2 * integralDiff;

    u = u1 + u2;
   

    oz=-2*a;
    y = Ct*q(:,k);


    dq(:,k) = [q(2,k);
               (g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(5,k));
               q(4,k);
               -q(9,k)*q(9,k)*sin(q(3,k))*cos(q(3,k))+(1/r)*(-((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(5,k)))*cos(q(3,k))-(-(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k)))*sin(q(8,k))*sin(q(3,k))+(((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k)))-g+g)*cos(q(8,k))*sin(q(3,k)));
               u(2,k)*sin(a)+u(1,k)*cos(a);
               q(7,k);
               -(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k));
               q(9,k);
               2*q(4,k)*q(9,k)*tan(q(3,k))+(1/r)*(1/cos(q(3,k)))*((-(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k)))*cos(q(8,k))+(((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k)))-g+g)*sin(q(8,k)));
               (u(2,k)*cos(a)-u(1,k)*sin(a))/cos(q(5,k));
               q(12,k);
               (g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k))-g;
                r0 - y];
    %αの更新
    da = -u(2,k)*cos(a)*tan(q(3,k))+u(1,k)*sin(a)*tan(q(3,k))+oz;

    % 1ステップ先の状態の更新
    q(:, k+1) = q(:,k) + dq(:,k)*Ts;
    a = a + da*Ts;
end

figure;
hold on;
plot(T,q(1,:),'DisplayName','x');
plot(T,q(6,:),'DisplayName','y');
plot(T,q(11,:),'DisplayName','z');
plot(T,q(3,:),'DisplayName','θ');
plot(T,q(8,:),'DisplayName','Φ');
legend('show');
xlabel('Time[s]');
ylabel('State [m][rad]');
yticks(min(ylim):0.1:max(ylim));
grid on;