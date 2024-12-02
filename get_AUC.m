function AUC=get_AUC(trace,center,Tau_front,Tau_back)

[preAmp, preAmpArg]=min(trace(center-Tau_front:center-1));
[postAmp, postAmpArg]=min(trace(center+1:center+Tau_back));

preAmpArg=preAmpArg+center-Tau_front-1;
postAmpArg=postAmpArg+center;

integral=sum(trace(preAmpArg:postAmpArg));
baseline=(preAmp+postAmp)/2*(postAmpArg-preAmpArg+1);

%AUC=(integral-baseline)/(postAmpArg-preAmpArg+1);
AUC=(integral-baseline);

end