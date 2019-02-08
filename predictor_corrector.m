% Analyzes positivity preservation of the Predictor-Corrector method
% applied to the Liouville-von Neumann equation.
%
% Copyright (c) 2018, Computational Photonics Group, Technical University 
% of Munich.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

clear;
close all;

% constants
e0 = 1.60217646e-19;
hbar = 1.05457168e-34;

% level count
N = 6;

% static electric field (V/m)
E = 9e9;

% Hamiltonian (diagonal elements)
H = zeros(N, N);
for n = 1:(N - 1)
    H(n + 1, n + 1) = H(n, n) + (1 - 0.1 * (n - 3)) * 2 * pi * 1e13 * hbar;
end

% dipole moment
d = 1e-29;

% Hamiltonian (off-diagonal elements)
for n = 1:(N-1)
    H(n, n + 1) = d * E;
    H(n + 1, n) = d * E;
end

% initial rho
rho_init = zeros(N, N);
rho_init(1, 1) = 1;

% enable higher precision
digits(32);
H = vpa(H);

% right-hand side Liouville-von Neumann equation
rhs = @(rho) ( -1i/hbar * (H * rho - rho * H));

% time span
te = 0.1e-12;
dt = 1e-16;
t = 0:dt:te;

% state variables
rho_me = rho_init;
rho_pc = rho_init;

% result: trace error
trace_me = zeros(size(t));
trace_pc = zeros(size(t));

% result: populations
pop_me = zeros(N, length(t));
pop_pc = zeros(N, length(t));

% prepare matrix exponential
U = expm(-1i * dt/hbar * H);

tic;
for i = 1:length(t)

    % Matrix exponential method
    rho_me = U * rho_me * U';

    % Predictor-Corrector method
    rho_pc = rho_pc + dt * rhs(rho_pc) ...
        + dt^2/2 * rhs(rhs(rho_pc)) ...
        + dt^3/4 * rhs(rhs(rhs(rho_pc))) ...
        + dt^4/8 * rhs(rhs(rhs(rhs(rho_pc))));

    trace_me(i) = trace(rho_me) - 1;
    pop_me(:, i) = diag(rho_me);

    trace_pc(i) = trace(rho_pc) - 1;
    pop_pc(:, i) = diag(rho_pc);
end
toc;

% display minimum population values
disp(['Minimum population values ME: ' num2str(min(min(real(pop_me))))]);
disp(['Minimum population values PC: ' num2str(min(min(real(pop_pc))))]);

% plot population rho_33
papersize = [ 15 12 ];
fig = figure('units', 'centimeters');
pos = get(gcf, 'pos');
set(gcf, 'pos', [pos(1) pos(2) papersize]);
plot(t/1e-12, real(pop_me(3, :)), '-.', 'Color', [0, 101, 189]/255,...
    'DisplayName', 'ME');
grid on;
hold on;
plot(t/1e-12, real(pop_pc(3, :)), '-', 'Color', [227, 114, 34]/255, ...
    'DisplayName', 'PC');
ax = gca;
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('Time/ps');
legend('show', 'Location', 'northeast');
ylabel('Population \rho_{33}/1');
ylim([-0.1 1]);

% create inset
axes(fig, 'Position', [0.25, 0.72, 0.4, 0.2]);
box on;
plot(t/1e-12, real(pop_me(3, :)), '-.', 'LineWidth', 1, 'Color', ...
    [0, 101, 189]/255, 'DisplayName', 'ME');
grid on;
hold on;
plot(t/1e-12, real(pop_pc(3, :)), '-', 'LineWidth', 1, 'Color', ...
    [227, 114, 34]/255, 'DisplayName', 'PC');
xlim([0.07 0.1]);
ylim([-0.1 0.6]);

set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize', papersize);
print(fig, 'pc.pdf', '-dpdf', '-fillpage');
