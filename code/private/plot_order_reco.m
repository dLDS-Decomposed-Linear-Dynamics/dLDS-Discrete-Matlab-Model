function plot_order_reco(x_vec)
figure;
x_vec = flip(x_vec);
num_x_vec_spec = 1;
for x_vec_spec = x_vec
x_vec_spec = x_vec_spec{1};
subplot(1,length(x_vec),num_x_vec_spec)
imagesc(x_vec_spec);
title(sprintf('Reconstruction from order %s',num2str(num_x_vec_spec)))
xlabel('Time');
num_x_vec_spec = num_x_vec_spec + 1;
end
end