function [N,M,q,Z] = getAlignmentFromFile(inputfile,stype)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & ...
        align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);
            for i=1:M
                counter = 0;
                for j=1:length(ind)
                    if( ind(j) )
                        counter = counter + 1;
                        Z(i,counter)=Molecule2Num( align_full(i).Sequence(j),stype );
                    end
                end
            end      
    q = max(max(Z));
end