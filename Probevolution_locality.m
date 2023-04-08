function [Trajectory_amino,Trajectory_nucleo, H, sitecount,sustcount,mutatedsites,Timeline,SustTimeline]=Probevolution_locality(S_amino,S_nucleo,e,h,T,M)
% % % % % % % % % % % % % % % % %
% S = the sequence that you start with to evolve.  Must be in number format
% e = family couplings
% h = local fields
% T = Selection temperature
% M = the number of evolutionary steps you want to run the simulation for.
%%%%%%%OUTPUTS%%%%%
% Trajectory = All the sequences produced by the simulation
% H = the hamiltonians for each sequence
% sitecount = the number of times each site got picked for potential mutation
% sustcount = the number of times each site actually got mutated
% Timeline = (MxN array) marks with a one the generation and site where a
% sampling was accepted. This can be the same aminoacid sequence
% SustTimeline = (MxN array) marks with a one the generation and site where
% a sustitution was accepted leading to a change in the aminoacid identity.
% mutated sites = If flags which sites chnaged during the simulation at
% least once.
    e=e/T;
    h=h/T:
    N=size(S_amino,2);
    gaps=S_amino==1;
    Trajectory_amino=zeros(M,N);
    Trajectory_nucleo=zeros(M,3*N);
    Trajectory_amino(1,:)=S_amino;
    Trajectory_nucleo(1,:)=S_nucleo;
    sitecount=zeros(1,N);
    sustcount=zeros(1,N);
    flag=zeros(M,1);
    Timeline=zeros(M,N);
    SustTimeline =zeros(M,N);
    for generation=2:M
        %Choose site
        is_gap=true;
        while is_gap % this loop continues checking if is_gap is true by comparing it to the site that was chosen.  If the chosen site is not a gap, it will be false and it will exit the loop. 
            ms=randvar(cumsum(ones(1,N)/N),1);
            is_gap=gaps(ms);
        end
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
            counter=counter+1; 
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
        Timeline(generation,ms)=1;
        Trajectory_amino(generation,:)=Sprueba;
        Trajectory_nucleo(generation,:)=Sprueba_nucleo;
        if (Trajectory_amino(generation,ms)~=Trajectory_amino(generation-1,ms))
            sustcount(ms)=sustcount(ms)+1;
            flag(generation)=1;
            SustTimeline(generation,ms)=1;
        end
    end

    mutatedsites=find(sustcount~=0);
    H=Generalhamiltonian(Trajectory_amino,e*T,h*T,2,1);
    end
    
    function X=randvar(P,n)
    % returns a realization for a discreate random variable given its
    % cumulative probability distribution.
    [~,X] = histc(rand(1,n),P);
    X=X+1;
    end
