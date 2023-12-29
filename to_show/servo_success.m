clear
format compact
close all

%x=[x,xd,θ,θd,β,y,yd,Φ,Φd,γ,z,zd]'
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

% 係数行列を拡大系にする.
At = [A ,zeros(12, 1); % 13*13
     -C             ];
Bt = [B; zeros(1, 3)];  % 13*3

%重みの決定
Q = diag([1, ...x
          1, ...xd
          1, ...θ
          100,...θd
          1, ...β
          1, ...y
          1,...yd
          1, ...Φ
          100,...Φd
          1, ...γ
          1, ...z
          1, ...zd
          1000]); %追従偏差 r0 - C*x
R=diag([1,1,1]);

%lqr法でゲインを計算
K = lqr(At, Bt, Q, R);

% シミュレーションの設定
Tf = 30;  % シミュレーション時間 [s]
Ts = 0.05; % サンプリング時間 [s]
T = 0:Ts:Tf;  % シミュレーションの時間軸
N = length(T);  % シミュレーションの時間ステップ数(=201)

% 初期状態の設定
x0 = [0, ...x
      0, ...xd
      0.1, ...θ
      0,...θd
      0, ...β
      0, ...y
      0,...yd
      -0.15, ...Φ
      0,...Φd
      0, ...γ
      0, ...z
      0, ...zd
      0]; %追従偏差 r0 - C*x

% 使用する変数等を0で宣言
x = zeros(13, N);
dx = zeros(13, N);
x(:, 1) = x0;
u = zeros(3, N); % 制御入力の初期化
x_ref = zeros(13, N);
u1 = zeros(3,1);
u2 = zeros(3,1);
a = 0;
da = 0;
y = zeros(13, N);

% 目標値
r0 = 1;

for k = 1:N-1

    % 状態フィードバックで制御入力を決定
    u(:,k) = -K*(x(:, k) - x_ref(:, k));

    oz=-2*a;

    y = r0 - C*x(:,k); % 追従偏差を取得

    dx(:,k) = [x(2,k);
               (g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*sin(x(5,k));
               x(4,k);
               -x(9,k)*x(9,k)*sin(x(3,k))*cos(x(3,k))+(1/r)*(-((g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*sin(x(5,k)))*cos(x(3,k))-(-(g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*sin(x(10,k))*cos(x(5,k)))*sin(x(8,k))*sin(x(3,k))+(((g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*cos(x(5,k))*cos(x(10,k)))-g+g)*cos(x(8,k))*sin(x(3,k)));
               u(2,k)*sin(a)+u(1,k)*cos(a);
               x(7,k);
               -(g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*sin(x(10,k))*cos(x(5,k));
               x(9,k);
               2*x(4,k)*x(9,k)*tan(x(3,k))+(1/r)*(1/cos(x(3,k)))*((-(g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*sin(x(10,k))*cos(x(5,k)))*cos(x(8,k))+(((g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*cos(x(5,k))*cos(x(10,k)))-g+g)*sin(x(8,k)));
               (u(2,k)*cos(a)-u(1,k)*sin(a))/cos(x(5,k));
               x(12,k);
               (g+0.5*g*(x(5,k)*x(5,k)+x(10,k)*x(10,k))+u(3,k))*cos(x(5,k))*cos(x(10,k))-g;
                y]; % 拡大系として追加した追従偏差の行

    

    %αの更新
    da = -u(2,k)*cos(a)*tan(x(3,k))+u(1,k)*sin(a)*tan(x(3,k))+oz;

    % 1ステップ先の状態の更新
    x(:, k+1) = x(:,k) + dx(:,k)*Ts;
    a = a + da*Ts;
end

figure;
hold on;
plot(T,x(1,:),'DisplayName','x');
plot(T,x(6,:),'DisplayName','y');
plot(T,x(11,:),'DisplayName','z');
plot(T,x(3,:),'DisplayName','θ');
plot(T,x(8,:),'DisplayName','Φ');
legend('show');
xlabel('Time[s]');
ylabel('State [m][rad]');
yticks(min(ylim):0.1:max(ylim));
grid on;