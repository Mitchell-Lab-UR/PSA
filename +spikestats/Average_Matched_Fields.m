function tinfo = Average_Matched_Fields(tinfo, minfo)
% function tinfo = Average_Matched_Fields(tinfo, minfo)
%
%  inputs:  tinfo - larger struct to append a new field
%           minfo - cell of M+1 structs with any set of fields
%                 - it will search all fields, average them into one
%                 - and then append field into tinfo
%           Note: ignore the 1st item of fields (will be empty, or
%                   drawn from an unmatched set)
%  outputs:  tinfo - now with appended struct

   M = length(minfo);
   for mm = 1:M
       mfields = fields(minfo);
       mfields
   end
   
end

