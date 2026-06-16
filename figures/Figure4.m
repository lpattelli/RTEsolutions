function Figure4
% FIGURE4 - CW anisotropic pole, ratio, and residue maps
%
% Description:
%   Generates panels 4a-c: the pole-matched parameter a, the corresponding
%   asymptotic RTE/DE ratio, and the incompatible pole- and residue-matched
%   diffusion coefficients in the anisotropic CW regime.
%
% Authors: Andre Liemert (ILM-ULM)
%          Lorenzo Pattelli (INRIM)
%          Fabrizio Martelli (UNIFI)
% License: MIT

  g_vals = linspace(0, 0.99, 34);
  eta_vals = logspace(-5, 0, 51);
  musp = 1;
  N = 101;

  [a_pole, ratio, g_vals, mua_vals, D_pole, D_res] = compute_a(g_vals, eta_vals, musp, N);
  logeta = log10(mua_vals / musp);
  [G, E] = meshgrid(g_vals, logeta);

  figure('Name', 'Figure 4: CW anisotropic asymptotic maps');

  %% Figure 4a - Pole-matched a
  subplot(1, 3, 1);
  draw_surface_with_floor_contours(G, E, a_pole.', [0.25, 0.35, 0.45, 0.55], 0);
  xlabel('g');
  ylabel('log_{10}(\mu_a/\mu_s'')');
  zlabel('a');
  title('Figure 4a');
  zlim([0, 0.8]);
  view([-63, 9]);

  %% Figure 4b - Asymptotic ratio at D_pole
  subplot(1, 3, 2);
  draw_surface_with_floor_contours(G, E, ratio.', [0.9, 0.99, 0.999, 0.9999], 0);
  xlabel('g');
  ylabel('log_{10}(\mu_a/\mu_s'')');
  zlabel('lim \phi/\psi');
  title('Figure 4b');
  zlim([0, 1]);
  view([-63, 9]);

  %% Figure 4c - Pole and residue branches
  subplot(1, 3, 3);
  hold on;
  surf(G, E, D_pole.', 'EdgeColor', 'none', 'FaceAlpha', 0.72);
  surf(G, E, D_res.', 'EdgeColor', 'none', 'FaceColor', [1.0, 0.55, 0.55], 'FaceAlpha', 0.35);
  add_floor_contours(G, E, D_pole.', [0.3, 0.33, 0.333, 0.3333], 0.2);
  xlabel('g');
  ylabel('log_{10}(\mu_a/\mu_s'')');
  zlabel('D\mu_s''');
  title('Figure 4c');
  zlim([0.2, 0.7]);
  view([-63, 9]);
  grid on;
  hold off;
end


function draw_surface_with_floor_contours(G, E, Z, levels, zfloor)
  surf(G, E, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.72);
  hold on;
  add_floor_contours(G, E, Z, levels, zfloor);
  grid on;
  hold off;
end


function add_floor_contours(G, E, Z, levels, zfloor)
  C = contourc(G(1,:), E(:,1), Z, levels);
  idx = 1;
  while idx < size(C, 2)
    npts = C(2, idx);
    xs = C(1, idx+1:idx+npts);
    ys = C(2, idx+1:idx+npts);
    plot3(xs, ys, zfloor*ones(1, npts), 'k-', 'LineWidth', 0.8);
    idx = idx + npts + 1;
  end
end
