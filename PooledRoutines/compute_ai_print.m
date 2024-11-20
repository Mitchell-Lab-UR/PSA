function [amed1,amed2] = compute_ai_print(ai1,ai2,topicname)

    amed1 = nanmedian(ai1);
    amed2 = nanmedian(ai2);
    med1 = (1+amed1)/(1-amed1);
    med2 = (1+amed2)/(1-amed2);
    med12 = mean([med1,med2]);
    [p1,h1] = signrank(ai1);
    [p2,h2] = signrank(ai2);
    [p12,h12] = ranksum(ai1,ai2);
    %******** disp stats in print
    disp([topicname,' summary statistics']);
    disp(sprintf('Marmo E: med(%7.4f) ai(%7.4f) p(%10.8f)',med1,amed1,p1));
    disp(sprintf('Marmo M: med(%7.4f) ai(%7.4f) p(%10.8f)',med2,amed2,p2));
    disp(sprintf('Marmo E vs M differ: med(%7.4f) p(%10.8f)',med12,p12));
    disp('raw');
    p1
    p2
    p12
    disp('  ');
    

return;
