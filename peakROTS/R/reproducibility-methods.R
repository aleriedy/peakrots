#----------------------------------
# OVERLAP METHODS
#----------------------------------

createUnion.new <- function(peaks){
	# Assumptions:
	# peaks is a matrix
	# peaks[,1] = chromosome name
	# peaks[,2] = peak start
	# peaks[,3] = peak end
	# no other columns
	
	# Next transform is useless if createUnion is only called from calculateUnionOverlap().
	peaks <- data.frame(chr = as.character(peaks[, 1]), start = as.numeric(peaks[, 2]), end = as.numeric(peaks[, 3]), stringsAsFactors=F);
	
	if (nrow(peaks) > 1) { # Has more than one row
	   
		cc <- sort(unique(peaks[, 1])); # Should be a short list.
		U <- NULL;
		
		## Evaluate each chromosome separately.
		for(i in 1:length(cc)) {
		
			## Fetch all peak regions in the particular chromosome and order them by start position.
			A <- peaks[peaks[,1] == cc[i], ];
			A <- A[order(A[, 2]), ];
			
			# Find the set of islands and add the to the union.
			if (nrow(A) >= 1) {
				a <- cummax(A[, 3]);
				start <- c(1, which(A[-1, 2] - a[-length(a)] > 0) + 1);
				U <- rbind(U, data.frame(
					chr   = rep(cc[i], length(start)),
					start = A[start, 2],
					end   = a[c(start[-1] - 1, length(a))], stringsAsFactors=F));
			}
		}
	} else U <- peaks;
	row.names(U) <- NULL;
	return(U);
}


calculateUnionOverlap <- function(peaks1, peaks2) {
	  
	peaks1 <- data.frame(chr = as.character(peaks1[, 1, drop = F]), start = as.numeric(peaks1[, 2, drop = F]), end = as.numeric(peaks1[, 3, drop = F]), stringsAsFactors=F);  
	peaks2 <- data.frame(chr = as.character(peaks2[, 1, drop = F]), start = as.numeric(peaks2[, 2, drop = F]), end = as.numeric(peaks2[, 3, drop = F]), stringsAsFactors=F);
	  
	U <- createUnion.new(rbind(peaks1, peaks2));
	  
	if (nrow(U) >= 1) {
		m1 <- rep(0, length = nrow(U));
		m2 <- rep(0, length = nrow(U));
		m <- rep(0, length = nrow(U));
		cc <- sort(unique(U[, 1]));
		
		## Index offset for chromosome groups
		ind <- 0;
		
		## Evaluate each chromosome separately.
		for (i in 1:length(cc)) {

			## All peak regions in the particular chromosome in the union list.
			UU <- U[U[, 1] == cc[i], ];

			for (j in 1:2) { # Same thing for both peaklists.
				
				## Catenate union and all peak regions in the particular chromosome in one of the lists.
				B <- {if (j == 1) peaks1[peaks1[, 1, drop = F] == cc[i], , drop = F] else peaks2[peaks2[, 1, drop = F] == cc[i], , drop = F]};
				A <- rbind(cbind(UU, index = 1:nrow(UU)), cbind(B, index = rep(NA, nrow(B))));

				## Order peaks by their start position.
				A <- A[order(A[, 2]), ];
				
				# Find the set of islands and add the to the union.
				if (nrow(A) >= 1) {
					a <- cummax(A[, 3]);
					start <- c(1, which(A[-1, 2] - a[-length(a)] > 0) + 1);
					
					# Remove one element cases. These don't contain items from individual lists.
					temp <- setdiff(start, c(start[-1] - 1, length(a)));
					
					# Flag found elements.
					m[ind + A[temp, 4]] <- ifelse(j == 1, 1, m[ind + A[temp, 4]] + 1);
					if (j == 1) m1[ind + A[temp, 4]] <- 1
					else m2[ind + A[temp, 4]] <- 1;
				}
			}
			ind <- ind + nrow(UU);
		}
	}
	
	## Calculate the proportion of peak regions in the union list that
	## share overlapping coordinates with both of the two separate peak lists.
	return(sum(m1 + m2 == 2) / length(m1));
}

 
