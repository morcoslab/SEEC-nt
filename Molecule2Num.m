function x=Molecule2Num(a,stype)
	switch (stype)
	case 1
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
	case 2
    switch(a)
        % full AA alphabet
        case 'A'
             x=1;
        case 'C'    
            x=2;    
        case 'G'    
            x=3;
        case 'T'
            x=4;
        case 'U'
            x=4;
        case '-'
            x=5;
        otherwise
            x=5;
    end
    end
end