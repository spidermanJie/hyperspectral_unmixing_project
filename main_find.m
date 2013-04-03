w_size = 10; 
time = 1:length(prices(1))
avg_price = filter(ones(w_size)./w_size,prices(stock_ind,:))
delta = avg_price(1:(end-1))-avg_price(2:end);
turning = find(delta(1:(end-1)).*delta(2:end)<0); 