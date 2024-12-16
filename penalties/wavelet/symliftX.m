function LS = symliftX(wname)
%SYMLIFT Symlets liftingX schemes.
%   LS = SYMLIFT(WNAME) returns the liftingX scheme specified
%   by WNAME. The valid values for WNAME are:
%      'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8'
%
%   A liftingX scheme LS is a N x 3 cell array such that:
%     for k = 1:N-1
%       | LS{k,1} is the liftingX "type" 'p' (primal) or 'd' (dual).
%       | LS{k,2} is the corresponding liftingX filter.
%       | LS{k,3} is the higher degree of the Laurent polynomial
%       |         corresponding to the previous filter LS{k,2}.
%     LS{N,1} is the primal normalization.
%     LS{N,2} is the dual normalization.
%     LS{N,3} is not used.
%
%   For more information about liftingX schemes type: lsinfoX.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Jun-2003.
%   Last Revision: 02-Nov-2007.
%   Copyright 1995-2007 The MathWorks, Inc.
%   $Revision: 1.1 $ $Date: 2012/02/08 09:52:45 $ 

Num = wstr2numX(wname(4:end));
switch Num
    %== sym2 ============================================================%
    case 2
        LS = {...
                'd',[-sqrt(3)],0; ...
                'p',[sqrt(3)-2 sqrt(3)]/4,1; ...
                'd',[1],-1 ...
            };
        LS(end+1,:) = {(sqrt(3)+1)/sqrt(2),(sqrt(3)-1)/sqrt(2),[]};
		% %--------------------  Num LS = 1  ----------------------% 
		% LS = {...                                                                 
		% 'd'                     [ -1.7320508075722079]                      [0]   
		% 'p'                     [ -0.0669872981075236  0.4330127018915160]  [1]   
		% 'd'                     [  0.9999999999994959]                      [-1]  
		% [  1.9318516525804916]  [  0.5176380902044105]                      []    
		% }; 
    case 2.1
		%--------------------  Num LS = 2  ----------------------% 
		LS = {...                                                                
		'd'                     [ -0.5773502691885155]                      [1]  
		'p'                     [  0.2009618943233436  0.4330127018926641]  [0]  
		'd'                     [ -0.3333333333327671]                      [0]  
		[  1.1153550716496254]  [  0.8965754721686846]                      []   
		};
    case 2.2        
		%--------------------  Num LS = 3  ----------------------% 
		LS = {...                                                                
		'd'                     [  0.5773502691900463]                      [0]  
		'p'                     [ -0.4330127018915159  2.7990381056783082]  [0]  
		'd'                     [ -0.3333333333332407]                      [1]  
		[  0.2988584907223872]  [  3.3460652149545598]                      []   
		};
    case 2.3                
		%--------------------  Num LS = 4  ----------------------% 
		LS = {...                                                                 
		'd'                     [  1.7320508075676158]                      [1]   
		'p'                     [ -0.4330127018926641 -0.9330127018941287]  [-1]  
		'd'                     [  0.9999999999980750]                      [2]   
		[ -0.5176380902041495]  [ -1.9318516525814655]                      []    
		};                                                                        
        
    %== sym3 ============================================================%
    case 3
		%--------------------  Num LS = 4  ----------------------% 
		LS = {...                                                                
		'd'                     [  0.4122865950085308]                      [0]  
		'p'                     [ -0.3523876576801823  1.5651362801993258]  [0]  
		'd'                     [ -0.4921518447467098 -0.0284590895616518]  [1]  
		'p'                     [  0.3896203901445617]                      [0]  
		[  0.5213212719156450]  [  1.9182029467652528]                      []   
		};                                                                              
    case 3.1
		%--------------------  Num LS = 2  ----------------------% 
		LS = {...                                                                
		'd'                     [ -0.4122865950517414]                      [1]  
		'p'                     [  0.4667569466389586  0.3523876576432496]  [0]  
		'd'                     [ -0.4921518449249469  0.0954294390155849]  [0]  
		'p'                     [ -0.1161930919191620]                      [1]  
		[  0.9546323126334674]  [  1.0475237290484967]                      []   
		}; 
    case 3.2
		%--------------------  Num LS = 5  ----------------------% 
		LS = {...                                                                 
		'd'                     [ -0.4122865950517414]                      [1]   
		'p'                     [ -1.5651362796324981  0.3523876576432496]  [0]   
		'd'                     [ -2.5381416988469603  0.4921518449249469]  [1]   
		'p'                     [  0.3896203899372190]                      [-1]  
		[  4.9232611941772104]  [  0.2031173973021602]                      []    
		};
    case 3.3
		%--------------------  Num LS = 6  ----------------------% 
		LS = {...                                                                 
		'd'                     [  2.4254972441665452]                      [1]   
		'p'                     [ -0.3523876576432495 -0.2660422349436360]  [-1]  
		'd'                     [  2.8953474539232271  0.1674258735039567]  [2]   
		'p'                     [ -0.0662277660392190]                      [-1]  
		[ -1.2644633083567955]  [ -0.7908493614571760]                      []    
		};                                                                        
        
    %== sym4 ============================================================%
    case 4
		%--------------------  Num LS = 16  ----------------------% 
		LS = {...                                                                
		'd'                     [  0.3911469419700402]                      [0]  
		'p'                     [ -0.1243902829333865 -0.3392439918649451]  [1]  
		'd'                     [ -1.4195148522334731  0.1620314520393038]  [0]  
		'p'                     [  0.4312834159749964  0.1459830772565225]  [0]  
		'd'                     [ -1.0492551980492930]                      [1]  
		[  1.5707000714496564]  [  0.6366587855802818]                      []   
		};                                                                       
    case 4.1
		%--------------------  Num LS = 19  ----------------------% 
		LS = {...                                                                 
		'd'                     [ -0.3911469419692201]                      [1]   
		'p'                     [  0.3392439918656564  0.1243902829339031]  [-1]  
		'd'                     [ -0.1620314520386309 -0.8991460629746448]  [2]   
		'p'                     [  0.4312834159764773 -0.2304688357916146]  [-1]  
		'd'                     [  0.6646169843776997]                      [2]   
		[  1.2500817546829417]  [  0.7999476804248136]                      []    
		};
 
    %== sym5 ============================================================%
    case 5
		%--------------------  Num LS = 29  ----------------------% 
		LS = {...                                                                
		'd'                     [ -0.9259329171294208]                      [0]  
		'p'                     [  0.4985231842281166  0.1319230270282341]  [0]  
		'd'                     [ -0.4293261204657586 -1.4521189244206130]  [1]  
		'p'                     [ -0.0948300395515551  0.2804023843755281]  [1]  
		'd'                     [  1.9589167118877153  0.7680659387165244]  [0]  
		'p'                     [ -0.1726400850543451]                      [0]  
		[  2.0348614718930915]  [  0.4914339446751972]                      []   
		};
    case 5.1
		%--------------------  Num LS = 23  ----------------------% 
		LS = {...                                                                
		'd'                     [  1.0799918455239754]                      [0]  
		'p'                     [  0.1131044403334987 -0.4985231842281165]  [1]  
		'd'                     [  2.4659476305614541 -0.5007584249312305]  [0]  
		'p'                     [ -0.0558424247659369 -0.2404034797205558]  [1]  
		'd'                     [  3.3265774193213002  1.3043080355478955]  [0]  
		'p'                     [ -0.1016623108755641]                      [0]  
		[ -2.6517078902691829]  [ -0.3771154446044534]                      []   
		};
    case 5.2
		%--------------------  Num LS = 24  ----------------------% 
		LS = {...                                                                
		'd'                     [  1.0799918455239754]                      [0]  
		'p'                     [  0.1131044403334987 -0.4985231842281165]  [1]  
		'd'                     [ -1.6937259364035369 -0.5007584249312305]  [0]  
		'p'                     [  0.2404034797205558 -0.0813027019760201]  [0]  
		'd'                     [  0.8958585825127051 -4.4713044896913088]  [1]  
		'p'                     [  0.1480132819787044]                      [0]  
		[  3.0742840094152357]  [  0.3252789907950670]                      []   
		};
        
    case 6
		%--------------------  Num LS = 1  ----------------------% 
		% Pow MAX = 0 - diff POW = 0
		%---+----+----+----+----+---%
		LS = {...                                                                 
		'd'                     [  0.2266091476053614]                      [0]   
		'p'                     [  1.2670686037583443 -0.2155407618197651]  [1]   
		'd'                     [ -0.5047757263881194  4.2551584226048398]  [-1]  
		'p'                     [ -0.0447459687134724 -0.2331599353469357]  [3]   
		'd'                     [ 18.3890008539693710 -6.6244572505007815]  [-3]  
		'p'                     [ -0.1443950619899142  0.0567684937266291]  [5]   
		'd'                     [  5.5119344180654508]                      [-5]  
		[ -1.6707087396895259]  [ -0.5985483742581210]                      []    
		};                                                                        

    case 7
		%--------------------  Num LS = 1  ----------------------% 
		% Pow MAX = 0 - diff POW = 0
		%---+----+----+----+----+---%
		LS = {...                                                                 
		'p'                     [  0.3905508237124110]                      [0]   
		'd'                     [ -0.3388639272262041  7.1808202373094066]  [0]   
		'p'                     [ -0.0139114610261505 -0.1372559452118446]  [2]   
		'd'                     [ 29.6887047769035310  0.1338899561610895]  [-2]  
		'p'                     [  0.1284625939282921 -0.0001068796412094]  [4]   
		'd'                     [ -7.4252008608107740 -2.3108058612546007]  [-4]  
		'p'                     [  0.0532700919298021  0.2886088139333021]  [6]   
		'd'                     [ -1.1987518309831993]                      [-6]  
		[  2.1423821239872392]  [  0.4667701381576485]                      []    
		};                                                                        

    case 8
		%--------------------  Num LS = 1  ----------------------% 
		% Pow MAX = 0 - diff POW = 0
		%---+----+----+----+----+---%
		LS = {...                                                                 
		'd'                     [  0.1602796165947262]                      [0]   
		'p'                     [  0.7102593464144563 -0.1562652322408773]  [1]   
		'd'                     [ -0.4881496179387070  1.8078532235524318]  [-1]  
		'p'                     [  1.7399180943774144 -0.4863315213006700]  [3]   
		'd'                     [ -0.5686365236759819 -0.2565755576271975]  [-3]  
		'p'                     [ -0.8355308510520870  3.7023086183759020]  [5]   
		'd'                     [  0.5881022226370752 -0.3717452749902822]  [-5]  
		'p'                     [ -2.1580699620177337  0.7491890598341392]  [7]   
		'd'                     [  0.3531271830147090]                      [-7]  
		[  0.4441986800900797]  [  2.2512448704197152]                      []    
		};                                                                        
        
    otherwise
        error('Wavelet:FunctionInput:Invalid_ArgVal',...
            'Invalid wavelet number.')
        
end