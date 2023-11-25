if Problem == 1
    K = 1;
    SIDE = "TOP";
    if SIDE == "TOP"
        A = 2;  B = 4;
    else
        A = 1;  B = 3;
    end

    if IEL == 2
        if SIDE == "TOP"
            STARTING = 11;      ENDING = 15;
        else
            STARTING = 1;       ENDING = 1;
        end
    elseif IEL == 1
        if SIDE == "TOP"
            STARTING = 51;      ENDING = 50;
        else
            STARTING = 1;       ENDING = 10;
        end
    end

    PRESS_TABLE = zeros((ENDING - STARTING + 1)*IEL,2);
    
    for I = STARTING:ENDING
        X = MATRIX{I,2};
        PRESS = MATRIX{I,7};

        if IEL == 2
            PRESS_TABLE(K,1) = X(A);
            PRESS_TABLE(K+1,1) = X(B);
            PRESS_TABLE(K,2) = PRESS(A);
            PRESS_TABLE(K+1,2) = PRESS(B);
            K = K+2;
        elseif IEL == 1
            PRESS_TABLE(K,1) = X;
            PRESS_TABLE(K,2) = PRESS;
            K = K+1;
        end
    end
    PRESS_TABLE
elseif Problem == 2
    if NY > 9
        K = 1;
        if IEL == 2
            STARTING = 73;      ENDING = 80;
        elseif IEL == 1
            STARTING = 305;     ENDING = 320;
        end

        PRESS_TABLE = zeros((ENDING-STARTING+1)*IEL,2);

        for I = STARTING:ENDING
            X = MATRIX{I-STARTING+1,2};
            PRESS = MATRIX{I-STARTING+1,7};

            if IEL == 2
                PRESS_TABLE(K,1) = X(2);
                PRESS_TABLE(K+1,1) = X(4);
                PRESS_TABLE(K,2) = PRESS(2);
                PRESS_TABLE(K+1,2) = PRESS(4);
                K = K+2;
            elseif IEL == 1
                PRESS_TABLE(K,1) = X;
                PRESS_TABLE(K,2) = PRESS;
            end
        end
        PRESS_TABLE
    else
        disp ('Pressure will be calculated only for 8*10-Q9 and 16*20-L4 case :)');
    end
end
