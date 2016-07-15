%% Dice Simulation
% Simulates 2 dices
N = 10000;
factor = zeros(1,N);
result = zeros(1,12);

for j = 1:N
  factor(j) = factorial(j);
  dice = randi(6);
  dice = dice + randi(6);
  result(dice) = result(dice) + 1;  
end

bar(result,'DisplayName', '2 Dice Simulation'); figure(gcf);

%% Choice change simulation
% Simulates choice
firstchoice = 0;
success = 0;
N = 1000;
T = 10000;

for j = 1:T    
    boxes = zeros(1,N); %initiate boxes
    prize = randi(N);   
    boxes(prize) = 1;   %set the prize in a random box

    choice = randi(N);  %choose the box

    if (boxes(choice) == 1) 
        firstchoice = firstchoice + 1;  %the first guess is correct
    end
end  % end of simulation without switching

for j = 1:T
    boxes = zeros(1,N); %initiate boxes
    prize = randi(N);   
    boxes(prize) = 1;   %set the prize in a random box

    choice = randi(N);  %choose the box
    
    if (choice ~= prize)  %if the box without prize is chosen, new boxes are 
                          % the one with the prize and the chosen one  
        random_choice = prize;                  
       
    else        % if the box with the prize is chosen
                %first, choose another random box that has no prize
        random_choice = randi(N); 
        while (random_choice == prize)
            random_choice = randi(N);
        end
       %now new boxes are the box with the prize and random box without prize         
    end
       boxes2 = [boxes(choice),boxes(random_choice)]; 
        %initial choice is first one
    
    
    if (boxes2(2) == 1)  %switch choice to the second
        success = success+1;
    end
    
end
f1sr = firstchoice/T*100;
f2sr = success/T *100;

f1sr
f2sr


%% Simple choice simulation
% Simulates simple choice

boxes = [0,1];
success = 0;

for j = 1:10000
    if (boxes(randi(2)) == 1) 
        success = success + 1;
    end
end

fssr = success /j *100;
fssr
