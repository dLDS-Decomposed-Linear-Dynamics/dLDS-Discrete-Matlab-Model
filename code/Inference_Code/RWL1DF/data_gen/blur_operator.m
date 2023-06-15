function frame_out = blur_operator(frame_in)

[M, N] = size(frame_in);
frame_out = zeros(M/2, N/2);

for kk = 1:M/2
    for nn = 1:N/2
        frame_out(kk, nn) = mean(reshape(frame_in(2*(kk-1)+1:2*kk, 2*(nn-1)+1:2*nn), [], 1));
    end
end


end