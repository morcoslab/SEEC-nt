function [isReplaceable,Replacing_NTseq]=isNB(CurrentNucleo,SubsAmino)
Replacing_NTseq=CurrentNucleo;
if SubsAmino == 1
    isReplaceable=false;
else
    amino_letters ='-ACDEFGHIKLMNPQRSTVWY';
    nucleo_letters ='ACGT-';



    NewResidue = amino_letters(SubsAmino);
    Codes4New=reshape(cell2mat(getfield(revgeneticcode,NewResidue)),3,[])';
    CurrentResidue = nucleo_letters(CurrentNucleo);

    Distances=3*squareform(pdist([double(CurrentResidue);double(Codes4New)],'hamming'));
    Distances=Distances(2:end,1);

    PossibleCodes4New=Codes4New(Distances<2,:);

    if isempty(PossibleCodes4New)
        isReplaceable=false;
    else
        isReplaceable=true;
        NewNTseq=PossibleCodes4New(randi(size(PossibleCodes4New,1)),:); %QUESTION: have we culled the possible codons for the new amino acid down to only those that are neighbors of the original codon?
        for i=1:3
            Replacing_NTseq(1,i)=Molecule2Num(NewNTseq(i),2);
        end
    end
end