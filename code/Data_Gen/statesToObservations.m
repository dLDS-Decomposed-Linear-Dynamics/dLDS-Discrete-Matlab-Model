function dFF = statesToObservations(D, latentX, sigVar)



dFF = cell(numel(latentX),1);
for ll = 1:numel(latentX)
    dFF{ll} = (D*latentX{ll}).' + sqrt(sigVar)*randn(size(latentX{ll},2), size(D,1));
    dFF{ll} = dFF{ll}.';
end

end