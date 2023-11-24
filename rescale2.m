function normtr=rescale2(tr,dim)
normtr=(tr-min(tr,[],dim))./(max(tr,[],dim)-min(tr,[],dim));
end