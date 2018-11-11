function [stepLength,newErr] = calculateStepLength(model,trueRec,source,dims,gradient,oldErr)
             
      

            
a(1)=1;
b(1)=251;
r=(sqrt(5)-1)/2;
d(1)=(b(1)-a(1))*(1-r);
alpha(1)=a(1)+d(1);
beta(1)=b(1)-d(1);
falpha=EvaluateSL(alpha(1),trueRec,gradient,model,source,dims);
fbeta=EvaluateSL(beta(1),trueRec,gradient,model,source,dims);

n=50;


for i = 2:n
 if falpha <= fbeta
    % New endpoints
   a(i) = a(i-1);                     
   b(i) = beta(i-1);                 
   % New internal points
   alpha(i) = a(i)+(1-r)*(b(i)-a(i)); 
   beta(i) = alpha(i-1);              
   % One new function evaluation
   fbeta = falpha;                    
   falpha = EvaluateSL(alpha(i),trueRec,gradient,model,source,dims);              
 else
   % New endpoints
   a(i) = alpha(i-1);                
   b(i) = b(i-1);                     
   % New internal points
   alpha(i) = beta(i-1);              
   beta(i) = b(i)-(1-r)*(b(i)-a(i));  
   % One new function evaluation
   falpha = fbeta;                    
   fbeta = EvaluateSL(beta(i),trueRec,gradient,model,source,dims);                
 end
 if (falpha || fbeta < oldErr) && (abs(alpha(i)-beta(i)) < 1)
     break 
 end

end
 %A=[alpha(i) beta(i)];
 %B=[falpha fbeta];
 stepLength=((alpha(i)+beta(i))/2);
 newErr=max((falpha+fbeta)/2);
end 
      
    
    
   
            
        
    
