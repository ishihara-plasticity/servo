% サーボ系線形状態方程式の導出
global As Bs1 Bs2 Cs;
delta = -a1*a2 + a3^2;
Cc = [1, 0, 0, 0]
Aa = [0, 0, 1, 0;
0, 0, 0, 1;
0, a3*a4/delta, a2*a5/delta, -a3*a6/delta;
0, -a1*a4/delta, -a3*a5/delta, a1*a6/delta]
Bb = [0, 0, -a2/delta, a3/delta]’
As = [ Aa, zeros(4,1);
-Cc, 0]
Bs1 = [Bb’, 0]’
Bs2 = [0, 0, 0, 0, 1]’
Cs = [Cc, 0]
Vs = ctrb(As, Bs1)
rank(Vs)
rank([A, B; Cc, 0])
% サーボ系設計
q1 = input(’q1: ’);
q2 = input(’q2: ’);
q3 = input(’q3: ’);
q4 = input(’q4: ’);
q5 = input(’q5: ’);
r = input(’r: ’);
Q = diag([q1, q2, q3, q4, q5])
R = [r]
global Fs;
[Fs, Ps, Es] = lqr(As, Bs1, Q, R);
Fs = -Fs
Es
% サーボ系シミュレーション
t0 = 0.0;
t1 = 20.0;
step = 0.001;
tspan = [t0:step:t1];
x0 = [0.1, 0.1, 0, 0]’;
eta0 = [0];
xs0 = [x0; eta0];
global ref;
ref = 0;
[Ts, Xs] = ode45(’nonlinservomodel’, tspan, xs0);
Us = Fs*Xs’;



figure(11);
plot(Ts, Xs);
legend(’x1’, ’x2’, ’x3’, ’x4’, ’eta’);
xlabel(’time [s]’);
ylabel(’state [rad,rad,rad/s,rad/s]’);
grid on;
figure(12);
plot(Ts, Us);
legend(’u’);
xlabel(’time [s]’);
ylabel(’input [Nm]’);
grid on;



function dxs = nonlinservomodel(t, xs)
global a1 a2 a3 a4 a5 a6
global Fs Cs
global ref
% 状態切出し
x = xs(1:4);
eta = xs(5);
y = Cs * xs;
% 重力加速度
g = 9.81;
% temporary 変数
x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
s = sin(x2);
c = cos(x2);
h = -a2*a2*s*s + a3*a3*c*c - a1*a2;
% 非線形状態方程式
fx =...
[x3;
x4;
(a2*a3*s*c^2*x3^2 - a2*a3*s*x4^2 + a2*a5*x3 - a3*a6*c*x4 +...
2*a2^2*s*c*x3*x4 + a3*a4*s*c)/h;
(-(a1+a2*s^2)*a2*s*c*x3^2 + a3^2*s*c*x4^2 - a3*a5*c*x3 +...
(a1+a2*s^2)*a6*x4 - 2*a2*a3*s*c^2*x3*x4 - (a1+a2*s^2)*a4*s)/h];
gx = [0;
0;
-a2/h;
a3*c/h];
us = Fs*xs;
dx = fx + gx*us;
% 目標値の変更
if t < 10
ref = 0.0;
else
ref = 1.0;
end
deta = ref - y;
dxs = [dx; deta];