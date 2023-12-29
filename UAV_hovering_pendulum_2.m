%q=[x,xd,θ,θd,β,y,yd,Φ,Φd,γ,z,zd]'

% パラメータ
%r = 0.6;% λ=2/3Lp
r = 2/3;
g = 9.8; % 重力加速度 [m/s^2]

A = [0 1 0 0 0 0 0 0 0;
    0 0 g 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 -g/r 0 0 0;
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0  0;];

Bl = [0 0 0;
    0 0 0;
    1 0 0;
    0 0 0;
    0 0 0;
    0 1 0;
    0 0 0;
    0 0 0;
    0 0 1;];

%重みの決定
Q = diag([1, 1, 1, 1, 1, 1, 1, 1, 1]);
R = diag([1,1,1]);

%lqr法
[K, S, E] = lqr(A, Bl, Q, R); % K がフィードバックゲイン行列

% % シミュレーションの設定
% Tf = 30; % シミュレーション時間 [s]
% Ts = 0.001; % サンプリング時間 [s]
% T = 0:Ts:Tf;  % シミュレーションの時間軸
% N = length(T);  % シミュレーションの時間ステップ数
% 
% % 目標値の定義
% %q_ref=[x,xd,θ,θd,β,y,yd,Φ,Φd,γ,z,zd]' 振子有りの場合
% q_ref = [0; 0; 0; 0; 0; 0; 1; 0];
% 
% % 初期状態
% q0 = [0; 0; 0; 0; 0; 0; 0; 0];
% 
% % 制御入力の初期化
% u = zeros(3, N);
% 
% % シミュレーションの実行
% q = zeros(8, N);
% dq = zeros(8, N);
% q(:, 1) = q0;
% for k = 1:N-1
% 
%     % 状態フィードバックによる制御入力の計算
%     u(:,k) = -K*(q(:,k) - q_ref);
% 
%     %状態方程式の計算
%     dq(:,k) = [q(2,k);
%                (g+0.5*g*(q(3,k)*q(3,k)+q(6,k)*q(6,k))+u(3,k))*sin(q(3,k));
%                u(1,k);
%                q(5,k);
%                -(g+0.5*g*(q(3,k)*q(3,k)+q(6,k)*q(6,k))+u(3,k))*sin(q(6,k))*cos(q(3,k));
%                u(2,k)/cos(q(3,k));
%                q(8,k);
%                (g+0.5*g*(q(3,k)*q(3,k)+q(6,k)*q(6,k))+u(3,k))*cos(q(3,k))*cos(q(6,k))-g];
% 
%     % 1ステップ先の状態の更新
%     q(:, k+1) = q(:,k) + dq(:,k)*Ts;
% end
% 
% %結果のプロット
% figure;
% title('state ');
% hold on;
% plot(T,q(1,:),'DisplayName','x', 'LineWidth', 3);
% plot(T,q(4,:),'DisplayName','y', 'LineWidth', 2);
% plot(T,q(7,:),'DisplayName','z', 'LineWidth', 2);
% plot(T,q(3,:),'DisplayName','β', 'LineWidth', 2);
% plot(T,q(6,:),'DisplayName','γ', 'LineWidth', 2);
% legend('show');
% xlabel('Time [s]');
% ylabel('State [m] [rad]');
% grid on;
% xticks(0:2:Tf);
% 
% 
% figure;
% title('input u');
% hold on;
% plot(T,u(1,:),'DisplayName','ωx', 'LineWidth', 2);
% plot(T,u(2,:),'DisplayName','ωy', 'LineWidth', 2);
% plot(T,u(3,:),'DisplayName','a', 'LineWidth', 2);
% legend('show');
% xlabel('Time [s]');
% ylabel('State [m/s^2] [rad]');
% grid on;
% xticks(0:2:Tf);
% 
% figure;
% title('dq');
% hold on;
% plot(T,dq(1,:),'DisplayName','x', 'LineWidth', 2);
% plot(T,dq(2,:),'DisplayName','xd', 'LineWidth', 2);
% plot(T,dq(3,:),'DisplayName','β', 'LineWidth', 2);
% plot(T,dq(4,:),'DisplayName','y', 'LineWidth', 2);
% plot(T,dq(5,:),'DisplayName','yd', 'LineWidth', 2);
% plot(T,dq(6,:),'DisplayName','γ', 'LineWidth', 2);
% plot(T,dq(7,:),'DisplayName','z', 'LineWidth', 2);
% plot(T,dq(8,:),'DisplayName','zd', 'LineWidth', 2);
% xlabel('Time [s]');
% ylabel('State [m/s^2] [rad]');
% legend('show');
% grid on;
% xticks(0:2:Tf);
% 
% 
% figure;
% sgtitle('dq'); % 全体のタイトル
% 
% rows = 4;
% cols = 2;
% 
% for i = 1:8
%     subplot(rows, cols, i);
%     plot(T, dq(i,:), 'DisplayName', ['dq(' num2str(i) ',:'], 'LineWidth', 2);
%     xlabel('Time [s]');
%     ylabel('State [m/s^2] [rad]');
%     legend('show');
%     grid on;
%     xticks(0:2:Tf);
% end
