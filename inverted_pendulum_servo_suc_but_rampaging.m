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

C = [1 0 0 0 0 0 0 0 0 0 0 0];

At = [A, zeros(12, 1);
      -C, 0];
Bt = [B; zeros(1, 3)];

%重みの決定
%        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
%        [x,xd, θ,θd, β, y,yd, Φ,Φd,  γ,  z, zd, w]
% Q = diag([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 100]);
% Q = diag([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6000]); %成功 時間がかかるけど、収束する。
Q = diag([1, 50, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 6000]);
R=diag([1,1,1]);

%lqr法
Kt = lqr(At, Bt, Q, R);
K1 = Kt(:, 1:12);
K2 = Kt(:, 13);
% K3 = [10; 0; 0];

% シミュレーションの設定
Tf = 70;  % シミュレーション時間 [s]
Ts = 0.05; % サンプリング時間 [s]
T = 0:Ts:Tf;  % シミュレーションの時間軸
N = length(T);  % シミュレーションの時間ステップ数(=201)

% 目標値の定義
q_ref = zeros(13, N);

% 初期状態
q0 = [0; 0; 0.1;  0; 0; 0; 0; -0.15; 0; 0; 0; 0; 0];
q = zeros(13, N);
dq = zeros(13, N);
q(:, 1) = q0;
u = zeros(3, N); % 制御入力の初期化

% 目標値 [m]
r0 = 0.5;


u1 = zeros(3,1);
u2 = zeros(3,1);

for k = 1:N-1

    % 状態フィードバックによる制御入力の計算
    u1(:,k) = -K1*(q(1:12, k) - q_ref(1:12, k));
    u2(:,k) = -K2*(q(13, k) - q_ref(13, k));

    u = u1 + u2;


    %状態方程式の計算 13状態
    dq(:,k) = [q(2,k);
        (g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(5,k));
        q(4,k);
        -q(9,k)*q(9,k)*sin(q(3,k))*cos(q(3,k))+(1/r)*(-((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(5,k)))*cos(q(3,k))-(-(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k)))*sin(q(8,k))*sin(q(3,k))+(((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k)))-g+g)*cos(q(8,k))*sin(q(3,k)));
        u(1,k);
        q(7,k);
        -(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k));
        q(9,k);
        2*q(4,k)*q(9,k)*tan(q(3,k))+(1/r)*(1/cos(q(3,k)))*((-(g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*sin(q(10,k))*cos(q(5,k)))*cos(q(8,k))+(((g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k)))-g+g)*sin(q(8,k)));
        u(2,k)/cos(q(5,k));
        q(12,k);
        (g+0.5*g*(q(5,k)*q(5,k)+q(10,k)*q(10,k))+u(3,k))*cos(q(5,k))*cos(q(10,k))-g;
        r0 - q(13, k)]; % 追加した13行目


    % 1ステップ先の状態の更新
    q(:, k+1) = q(:,k) + dq(:,k)*Ts;
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
grid on;

% 入力をグラフ化
% figure;
% T2 = T(1:end-1);
% hold on;
% plot(T2, u(1,:), 'DisplayName', 'u1');
% plot(T2, u(2,:), 'DisplayName', 'u2');
% plot(T2, u(3,:), 'DisplayName', 'u3');
% title('Control Inputs Over Time');
% xlabel('Time [s]');
% ylabel('Control Input Value');
% legend('show');
% grid on;