function samp_out = samp_DCT(image_mat, M_DCT)

K = 0.5*(-1 + sqrt(1 + 8*M_DCT));

image_DCT = dct2(image_mat);
image_DCT = fliplr(image_DCT(1:K, 1:K));

samp_out = zeros(M_DCT, 1);

samp_out(1) = diag(image_DCT, K-1);
for kk = 2:K
    samp_out(sum(0:kk-1)+1:sum(0:kk-1)+kk) = diag(image_DCT, K-kk);
end


end
