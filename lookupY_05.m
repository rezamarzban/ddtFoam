% Windows, Cantera

clear all; %clc;
thisfilename = [mfilename('fullpath') '.m'];
write = 1;
showplots = 0;
dir = 'Z:\OpenFOAM\ettner-1.7.1\run\Cantera\';
writetofilename = 'cTable_fpT_07';
logfilename = [dir writetofilename '.log']; % full logfile name

fn = [dir writetofilename '.csv']; %full filename

h2=[0:1:30 32:2:60]*1e-2; nh2=length(h2); % H2 mole fraction 
p=1e5*[0.1 0.9 1.1 2 5 10:10:40 60 80 100 125 150]; np =length(p);
T=[250:20:1600]; nT=length(T);

fH = 0*h2;

len = nh2*np*nT;
disp(len);
reply = input('Continue? Y/N [Y]: ', 's');
if (reply=='N')
    error('program stopped by user')
end

mech = 'OConaire.cti';
gas0 = importPhase(mech); 
nsp = nSpecies(gas0);
ynames=speciesNames(gas0);

% find Hydrogen nitrogen, and oxygen indices
ih2 = speciesIndex(gas0,'H2');
io2  = speciesIndex(gas0,'O2');
in2  = speciesIndex(gas0,'N2');
ih2o = speciesIndex(gas0,'H2O');
ioh  = speciesIndex(gas0,'OH');
ih  = speciesIndex(gas0,'H');

x = zeros(nsp,1);
Yout = NaN(nsp,nh2,np,nT);

fig_num = 0; % 0=no plots

% log file:
if(write)
       logid = fopen(logfilename, 'w');
       %fprintf(logid, 'Induction times for constant-volume explosion\n');
       fprintf(logid, 'burned gas composition\n');
       fprintf(logid, 'fH, p, Tu');
       for(i=1:length(ynames))
         fprintf(logid,', ');
         fprintf(logid, ynames{i});
       end
       fprintf(logid,'\n');
       
end

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

               equilibrate(gas0,'HP')
               mtemp=massFractions(gas0);
               Yout(:,i,j,k)=mtemp;
               Tout(i,j,k) = temperature(gas0);
                
               
               if(write)
                  fprintf(logid, '%s, %s, %s %s\n', num2str(fH(i),'%1.5f'),num2str(p(j),'%1.3e'), ...
                      num2str(T(k),'%4.1f'),num2str(mtemp',', %1.5e'));
               end
            
        end
    end
end



if(write)
       fclose(logid); % close log file
disp('writing file');
%disp(fn);

%d = date;
	
fid = fopen(fn, 'w');
%fprintf(fid, 'Shock and Detonation Toolbox\n');
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


fprintf(fid, '\n\nfields\n');
fprintf(fid, '3\n');
fprintf(fid, '(\n');

fprintf(fid, '\n');
fprintf(fid, ' {\n');
fprintf(fid, '     name            fH;\n');
fprintf(fid, '     min             %1.5f;\n',min(fH));
fprintf(fid, '     max             %1.5f;\n',max(fH));
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
fprintf(fid, '     name            Tu;\n');
fprintf(fid, '     min             %4.1f;\n',min(T));
fprintf(fid, '     max             %4.1f;\n',max(T));
fprintf(fid, '     N               %d;\n',length(T)-1);
fprintf(fid, ' }\n');
%fprintf(fid, '\n');
fprintf(fid, ');\n\n');


fprintf(fid, '\noutput\n');
fprintf(fid, '%d\n',nsp);
fprintf(fid, '(\n');
for i=1:nsp

    fprintf(fid, ' {    name \t %s; }\n',ynames{i});
    % fprintf(fid, ' {    name \t O2; }');
end
fprintf(fid, ')\n');
fprintf(fid, ';\n\n');

fprintf(fid, 'values\n');        
fprintf(fid, '%d\n',3+nsp);
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

fprintf(fid, '\n%d  // Tu\n( ',len);
for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
    fprintf(fid, '\n%4.1f ',T(k));
        end
    end
end
fprintf(fid, '\n)\n');


for(spec=1:nsp)
  fprintf(fid, '\n%d //%s\n( ',len,ynames{spec});
for(i=1:nh2)    
    for(j=1:np)        
        for(k=1:nT)
            fprintf(fid, '\n%1.4e ',Yout(spec,i,j,k));
        end
    end
end
  fprintf(fid, '\n)\n');
end

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

%titletext=[mech(1:end-4) ', p = ' num2str(p0/1e5,'%1.2f') ' bar, T_0 = ' num2str(T0,'%1.0f') ' K'];
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
