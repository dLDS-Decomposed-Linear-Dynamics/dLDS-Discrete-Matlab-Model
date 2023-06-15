function image_mat = samp_DCT_trans(samp_out, M_DCT, image_size)

K = 0.5*(-1 + sqrt(1+8*M_DCT));
image_TMP = zeros(K); 

image_TMP(1, end) = samp_out(1);
for kk = 2:K 
    TMP = zeros(K);
    TMP(1:kk, end-kk+1:end) = eye(kk);
    image_TMP(TMP==1) = samp_out(sum(0:kk-1)+1:sum(0:kk-1)+kk);
end

image_DCT = zeros(image_size);
image_DCT(1:K, 1:K) = fliplr(image_TMP);

image_mat = idct2(image_DCT);
end

