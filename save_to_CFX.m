function save_to_CFX(geo)

    fp=fopen("jk_A.rtzt", "w"); % Mentés Meanline Data File Format-ba

    w=3; % lapátvastagság
    LEmeshExt = 0.5; % Háló hosszabításának aránya radiális irányban a belépő élnél
    TEmeshExt = 0.2; % Háló hosszabításának aránya radiális irányban a kilépő élnél
    slices = 5; % rétegek száma b2 hosszon
    
    fprintf(fp,"%d 0",geo.N_lapat); % Lapátszám; "number of splitters"
    fprintf(fp,"\n");
    fprintf(fp,"\n 0.00 %d N M", slices); % "pitch fraction"; rétegek száma b2 hosszon; N-lapátvastagság normális irányba; M-?
    
    LEr2 = sqrt(geo.x_g(1,1)^2+geo.y_g(1,1)^2); % Belépő él sugara
    LEr1 = LEr2 * LEmeshExt; % Háló legkisebb sugara

    TEr1 = sqrt(geo.x_g(end,1)^2+geo.y_g(end,1)^2); % Kilépő él sugara
    TEr2 = TEr1 * (1+TEmeshExt); % Háló legnagyobb sugara
    
    radNum = 5; % Hálóhosszabbítás felbontása
    LErads = linspace(LEr1,LEr2,radNum+1);
    TErads = linspace(TEr1,TEr2,radNum+1);

    LEtheta = acos(geo.x_g(1,1)/LEr2); % Theta szög a belépő élnél
    TEtheta = acos(geo.x_g(end,1)/TEr1); % Theta szög a kilépő élnél

    for i = 1:slices
        % Réteg magassága; "z" koordináta
        b = -1000 * geo.b2 / (slices-1) * (i-1);
        
        % "span fraction"-rétegek felosztása 0.0-tól 1.0-ig?;
        % ;rétegben megadott pontog száma
        fprintf( fp,"\n \n %.5f %d", 1.0/(slices-1) * (i-1), 2*radNum+length(geo.x_g(:,1)) );
        
        % Háló hosszabbított szakasza (lapátvastagság=0.0)
        for j = 1:radNum
            % r, theta, z, b
            fprintf(fp,"\n \t \t %.5f %.5f %.5f 0.00000", 1000*LErads(j), LEtheta, b);
        end
        
        % Lapát szakasza
        for j = 1:length(geo.x_g(:,1))
            r = sqrt(geo.x_g(j,1)^2+geo.y_g(j,1)^2);
            theta = acos(geo.x_g(j,1)/r);

            % r, theta, z, b
            fprintf(fp,"\n \t \t %.5f %.5f %.5f %.5f", 1000*r, theta, b, w);
            plot(r, theta, '*'), hold on
        end
        
        % Háló hosszabbított szakasza (lapátvastagság=0.0)
        for j = 1:radNum
            % r, theta, z, b
            fprintf(fp, "\n \t \t %.5f %.5f %.5f 0.00000", 1000*TErads(j+1), TEtheta, b);
        end

    end
    
    fclose(fp);
end
