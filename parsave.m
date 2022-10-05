function parsave(fname,epoch_struct_noHm,epoch_struct_Hm)
  save(fname,'-v7.3','epoch_struct_noHm', 'epoch_struct_Hm')
end