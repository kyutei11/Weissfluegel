function [target_alpha itl_fzero] = get_alphaTarget(targetCl, alpha0, com)
  com.targetCl = targetCl;
  [target_alpha itl_fzero]= fzero('get_dCl', alpha0, com, 1.e-6);
endfunction
