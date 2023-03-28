function [Trajectory_amino,Trajectory_nucleo, H, sitecount,sustcount,Timeline,Sustitutiontime,R,Sustitutiontimeline,mutatedsites]=Probevolution_locality(S_amino,S_nucleo,e,h,M)
% % % % % % % % % % % % % % % % %
% S = the sequence that you start with to evolve.  Must be in number format
% e = family couplings
% h = local fields
% M = the number of evolutionary steps you want to run the simulation for.
%%%%%%%OUTPUTS%%%%%
% Trajectory = All the sequences produced by the simulation
% H = the hamiltonians for each sequence
% sitecount = the number of times each site got picked for potential mutation
% sustcount = the number of times each site actually got mutated
% Timeline = [not entirely sure but I think this is the time each mutation
%       happened assuming a poissonian process. -CMN]
% Sustitutiontime =
% R = ?
% Sustitutiontimeline = same as timeline but shows the position where the
%       substitution happened
% mutated sites = ?
%to plot the hamitonian trajectory: figure();H2=plot(H)
N=size(S_amino,2);
gaps=S_amino==1;
Trajectory_amino=zeros(M,N);
Trajectory_nucleo=zeros(M,3*N);
Trajectory_amino(1,:)=S_amino;
Trajectory_nucleo(1,:)=S_nucleo;
sitecount=zeros(1,N);
sustcount=zeros(1,N);
flag=zeros(M,1);
Timeline=1:M;
Timeline=cumsum(Timeline);
coeff=M/50;
R=zeros(coeff,1);
Sustitutiontimeline=zeros(M,N);
for generation=2:M
    %Choose site
    is_gap=true;
    while is_gap % this loop continues checking if is_gap is true by comparing it to the site that was chosen.  If the chosen site is not a gap, it will be false and it will exit the loop. 
        ms=randvar(cumsum(ones(1,N)/N),1);
        is_gap=gaps(ms);
    end
    %ms_nucleotide=idivide(int16(ms),3)+1; %DOUBLE CHECK THIS! This seems to be indexing an amino acid position to 1/3 of a nucleotide position instead of 3x.
    ms_nucleotide=3*int16(ms)-2;
    isReplaceable=false; 
    counter=0;
    while not(isReplaceable)
        distP=siteprobdistribution(Trajectory_amino(generation-1,:),h,e,ms,2);
        alpha=int8(randvar(distP,1));
        if alpha == 1 
            continue 
        end
        [isReplaceable,Replacing_NTseq]=isNB(Trajectory_nucleo(generation-1,ms_nucleotide:ms_nucleotide+2),alpha);
        counter=counter+1; %QESTION: how is the counter being applied if you use a continue statment above this location. 
        if counter>100 %if you sample the prob distribution more than 100 times at this site and keep picking a gap, then you must go back to the top and pick another mutation site.
            isReplaceable=true;
            alpha=Trajectory_amino(generation-1,ms);
            Replacing_NTseq = Trajectory_nucleo(generation-1,ms_nucleotide:ms_nucleotide+2);
        end
    end
    Sprueba_nucleo=Trajectory_nucleo(generation-1,:);
    Sprueba_nucleo(ms_nucleotide:ms_nucleotide+2)=Replacing_NTseq; %update the nucleotide sequence with the whole new codon from the amino acid mutation
    Sprueba=Trajectory_amino(generation-1,:); %extract the new sequence from the trajectory
    Sprueba(ms)=alpha;
    sitecount(ms)=sitecount(ms)+1;
    Trajectory_amino(generation,:)=Sprueba;
    Trajectory_nucleo(generation,:)=Sprueba_nucleo;
    if (Trajectory_amino(generation,ms)~=Trajectory_amino(generation-1,ms))
        sustcount(ms)=sustcount(ms)+1;
        flag(generation)=1;
        Sustitutiontimeline(generation,ms)=Timeline(generation);
    end
    if (mod(generation,50)==0)
        Sustitutiontime=Timeline(flag==1);
        differ=diff(Sustitutiontime);
        R(generation/50)=var(differ)/mean(differ);
    end
end
Sustitutiontime=Timeline(flag==1);
mutatedsites=find(sustcount~=0);
H=Generalhamiltonian(Trajectory_amino,e,h,2,1);
end

function X=randvar(P,n)
% returns a realization for a discreate random variable given its
% cumulative probability distribution.
[~,X] = histc(rand(1,n),P);
X=X+1;
end
