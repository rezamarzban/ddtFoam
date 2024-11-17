% Windows, Cantera
% tabulation of log10(tIgn)

clear all; 
thisfilename = [mfilename('fullpath') '.m'];
write = 1;
showplots = 0;
dir = 'Z:\OpenFOAM\ettner-2.1.1\run\Cantera\';
writetofilename = 'LOG_tignTable_fpT_09';
logfilename = [dir writetofilename '.log']; % full logfile name

fn = [dir writetofilename '.csv']; %full filename

h2=[1 5 10 20 30 40 50 60]*1e-2; 
nh2=length(h2); % H2 mole fraction 
  fH = 0*h2;
  
p=1e5*[0.1 0.9 1.1 2 5 10 20 30 50 75 100]; np =length(p);


T=[800:20:2000]; nT=length(T);
len = nh2*np*nT;
disp(len);

mech = 'OConaire.cti';
gas0 = importPhase(mech); 
nsp = nSpecies(gas0);

% find Hydrogen nitrogen, and oxygen indices
ih2 = speciesIndex(gas0,'H2');
io2  = speciesIndex(gas0,'O2');
in2  = speciesIndex(gas0,'N2');
ih2o = speciesIndex(gas0,'H2O');
ioh  = speciesIndex(gas0,'OH');
ih  = speciesIndex(gas0,'H');

x = zeros(nsp,1);
tIgn=NaN(nh2,np,nT);

fig_num = 0; % 0=no plots

% log file:
if(write)
       logid = fopen(logfilename, 'w');
       fprintf(logid, 'Induction times\n');
       fprintf(logid, 'fH, p, T, tign, log10(tign)\n');
       
end

dtdT = 0; % maximum error dt_ign/dT

for(i=1:nh2)
    for(j=1:np)
        for(k=1:nT)
               disp([i j k]);
               x=zeros(nsp,1);
               x(ih2)=h2(i);
               o2(i)=0.21*(1-h2(i)); x(io2)=o2(i);
               n2(i)=1-h2(i)-o2(i);  x(in2)=n2(i);               
               set(gas0,'Temperature',T(k),'Pressure',p(j),'MoleFractions',x);
               mtemp=massFractions(gas0);
                fH(i)=mtemp(ih2);

               [out] = explosion(gas0,fig_num); % constant volume explosion
               %[out] = explosionHP(gas0,fig_num); % constant pressure explosion
               
               tIgn(i,j,k) = out.ind_time;    
               if(tIgn(i,j,k)>0.9) 
                   tIgn(i,j,k)=1e1;
               end  
               
               if(write)
                  fprintf(logid, '%s, %s, %s, %s, %s\n', num2str(fH(i),'%1.5f'),num2str(p(j),'%1.3e'), ...
                      num2str(T(k),'%4.1f'),num2str(tIgn(i,j,k),'%1.5e'),num2str(log10(tIgn(i,j,k)),'%1.5f'));
               end
               
               if(tIgn<0.1)
               if(k>1)
                  dtdT_temp = (tIgn(i,j,k) - tIgn(i,j,k-1))/(T(k)-T(k-1));
                  if(dtdT_temp>dtdT)
                      dtdT = dtdT_temp
                      dtdT_info = [i,j,k]
                  end
               end
               end
            
        end
    end
end



if(write)
       fclose(logid); % close log file
disp('writing file');
disp(fn);

	
fid = fopen(fn, 'w');
fprintf(fid, '/*--------------------------------*- C++ -*---------------------------------*\\\n');
fprintf(fid, '| =========                |                                                 |\n');
fprintf(fid, '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
fprintf(fid, '|  \\    /   O peration     | Version:  1.7.1                                 |\n');
fprintf(fid, '|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |\n');
fprintf(fid, '|    \\/     M anipulation  |                                                 |\n');
fprintf(fid, '\\*--------------------------------------------------------------------------*/\n');

fprintf(fid, '\nFoamFile\n');
fprintf(fid, '{\n');
fprintf(fid, 'version     2.0;\n');
fprintf(fid, 'format      ascii;\n');
fprintf(fid, 'class       dictionary;\n');
fprintf(fid, 'location    "constant";\n');
fprintf(fid, 'object      %s;\n',writetofilename);
fprintf(fid, '}\n');
fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
fprintf(fid, '\n// created with %s',thisfilename);
fprintf(fid, '\n// %s mechanism ',mech(1:end-4));


fprintf(fid, '\n\noptions\n');
fprintf(fid, '{\n');
fprintf(fid, '    logarithmic true;\n');
fprintf(fid, '}\n');

fprintf(fid, '\n\nfields\n');
fprintf(fid, '3\n');
fprintf(fid, '(\n');

fprintf(fid, '\n');
fprintf(fid, ' {\n');
fprintf(fid, '     name            fH;\n');
fprintf(fid, '     min             %1.5f; // xH2=%1.3f\n',min(fH),min(h2));
fprintf(fid, '     max             %1.5f; // xH2=%1.3f\n',max(fH),max(h2));
fprintf(fid, '     N               %d;\n',length(fH)-1);
fprintf(fid, ' }\n');
%fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, ' {\n');
fprintf(fid, '     name            p;\n');
fprintf(fid, '     min             %1.3e;\n',min(p));
fprintf(fid, '     max             %1.3e;\n',max(p));
fprintf(fid, '     N               %d;\n',length(p)-1);
fprintf(fid, ' }\n');
%fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, ' {\n');
fprintf(fid, '     name            T;\n');
fprintf(fid, '     min             %4.1f;\n',min(T));
fprintf(fid, '     max             %4.1f;\n',max(T));
fprintf(fid, '     N               %d;\n',length(T)-1);
fprintf(fid, ' }\n');
%fprintf(fid, '\n');
fprintf(fid, ');\n\n');


fprintf(fid, '\noutput\n');
fprintf(fid, '%d\n',1);
fprintf(fid, '(\n');
fprintf(fid, ' {    name \t tIgn; }');
fprintf(fid, '\n)\n');
fprintf(fid, ';\n\n');

fprintf(fid, 'values\n');        
fprintf(fid, '4\n');
fprintf(fid, '(\n');

fprintf(fid, '\n%d  // fH\n( ',len);

for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
    fprintf(fid, '\n%1.5f ',fH(i));
        end
    end
end
fprintf(fid, '\n)\n');

fprintf(fid, '\n%d  // p\n( ',len);
for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
    fprintf(fid, '\n%1.3e ',p(j));
        end
    end
end
fprintf(fid, '\n)\n');

fprintf(fid, '\n%d  // T\n( ',len);
for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
    fprintf(fid, '\n%4.1f ',T(k));
        end
    end
end
fprintf(fid, '\n)\n');


  fprintf(fid, '\n%d // log10(tIgn) \n( ',len);
for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
            fprintf(fid, '\n%1.5f ',log10(tIgn(i,j,k)));
        end
    end
end
  fprintf(fid, '\n)\n');

fprintf(fid, '\n)\n');
fprintf(fid, ';\n');

end  % write


if(showplots)
fs=20;
ms=8;
lw=3;
marker = 'd';

colori = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1; 0.8 0.8 0.8; 0.7 0 0; 0 0 0.7; 0 0.7 0]*0.8;
lines={'-','--','-.','.'};
figure1 = figure('Units','normalized','Position',[0.0 0.2 0.8 0.6],'PaperUnits','centimeters','PaperSize',[21 29.7],...
    'PaperOrientation','landscape','PaperType','A4','Color','w','PaperPosition',[2 2 16 10]);
axes1=axes('Parent',figure1,'FontSize',fs,'XScale','linear','YScale','log','Color','w');

for(i=1:nh2)
   
for(k=1:np)    
plot(T,tIgn(i,:,k),'DisplayName',['p_0 = ' num2str(p(k)/1e5,'%2.1f') ' bar'],'LineWidth',lw,'LineStyle',lines{i},'Marker','none','MarkerSize',ms,...
     'Color',colori(k,:));
hold on;
end

end

set(gca,'YScale','log');
yext=[1e-6 1e-0]; ylim(yext);
xext=[T(1) T(end)]; xlim(xext);
leg1=legend('show','Location','NorthEast');
set(leg1,'FontSize',fs-6);


titletext=['Ignition delay time (mechanism:  ' mech(1:end-4) ')'];
title(titletext);
ylabel('Ignition time t_i_g_n   [s]');
xlabel('Initial temperature T_0   [-]');

 h2=text(xext(2),yext(1),thisfilename,'VerticalAlignment','bottom','HorizontalAlignment','right',...
     'BackgroundColor',[1 1 1],'Interpreter','none');


hold off;

%   surf(T,h2 ,tIgn(:,:,1))
%   view(axes1,[136.5 42]); 
end
