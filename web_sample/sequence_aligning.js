var numArgChoice = 5; // A, T(U), C, G, invalid nucleotide 
var complimentaryMap = {'A':'T','C':'G','G':'C','T':'A','M':'K','K':'M','Y':'R','R':'Y','V':'B','B':'V','H':'D','D':'H',
						'a':'t','c':'g','g':'c','t':'a','m':'k','k':'m','y':'r','r':'y','v':'b','b':'v','h':'d','d':'h'};
/**
 * @param seq    sequence
 * @return an integer array
 */
function str_seq_to_num_code(seq){
	var codedSeq = [];
	codedSeq.length = seq.length;
	for (var i = 0; i < seq.length; ++i){
		switch(seq.charAt(i).toUpperCase()) {
			case 'A':
				codedSeq[i] = 0;
				break;
			case 'C':
				codedSeq[i] = 1;
				break;
			case 'G':
				codedSeq[i] = 2;
				break;
			case 'T':
				codedSeq[i] = 3;
				break;
			case 'U':
				codedSeq[i] = 3;
				break;
			default:
				codedSeq[i] = 4;
		}
	}
	return codedSeq;
}

function reverse_DNA_direction(original) {
	var after = '';
	for (var i = 0; i < original.length; i++){
		var aft = 'N';
		var ori_rev = original.charAt(original.length - i -1);
		aft = complimentaryMap[ori_rev];
		after += aft;
	}
	return after;
}

/**
 * Generate scoring matrix from match/mismatch score
 *
 * @param matchMismatchArray
 *
 * @return sqaure scoring matrix. The last row and column are zero, for
 * matching an ambiguous residue.
 */
function get_score_matrix(matchMismatchArray){
	matchScore = matchMismatchArray[0];
	mismatchScore = matchMismatchArray[1]
	var m = [];
	if (mismatchScore > 0) {
		console.warn("the mismatch score should be negative");
		mismatchScore = -mismatchScore; 
	} 
	for (var i = 0; i < numArgChoice - 1; ++i) {
		m[i] = [];
		m[i][numArgChoice - 1] = 0;
		for (var j = 0; j < numArgChoice - 1; ++j)
			m[i][j] = i == j? matchScore : mismatchScore;
	}
	m[numArgChoice-1] = [];
	for (var j = 0; j < numArgChoice; j++) m[numArgChoice-1][j] = 0;
	return m;
}

function align(target, query, matchMismatchArray, gapsc){
	// convert string bases to coded integers
	var t = str_seq_to_num_code(target);
	var query = str_seq_to_num_code(query);
	var matrix = get_score_matrix(matchMismatchArray);
	var qp = [];
	for (var j = 0; j < matrix.length; j++) {
		qp[j] = [];
		for (var i = 0; i < query.length; i++)
			qp[j][i] = matrix[j][query[i]];
	}

	// these are penalties which should be non-negative
	var gapo = Math.abs(Number(gapsc[0]));
	var gape = Math.abs(Number(gapsc[1]));
	var gapoe = gapo + gape;

	// set up
	var H = [], E = [], z = [], score, max = 0, end_i = -1, end_j = -1;
	H[0] = 0; 
	E[0] = -2*gapoe;
	for (var j = 1; j <= query.length; j++) {
		H[j] = H[0] - gape*(j - 1);
		E[j] = E[0] - gape*j;
	}
	// the DP loop
	for (var i = 0; i < t.length; ++i) {
		var m = 0, mj = -1;
		var qpi = qp[t[i]];
		z[i] = [];

		var h1 = -(gapoe + gape * i);
		var f  = -(2*gapoe + gape * i);

		for (var j = 0; j < query.length; j++) {
			// At the beginning of the loop: h=H[j]=H(i-1,j-1), e=E[j]=E(i,j), f=F(i,j) and h1=H(i,j-1)
			// If we only want to compute the max score, delete all lines involving direction "d".
			var e = E[j], h = H[j], d;
			H[j] = h1;           // set H(i,j-1) for the next row
			h += qpi[j];         // h = H(i-1,j-1) + S(i,j)
			d = h >= e? 0 : 1;
			h = h >= e? h : e;
			d = h >= f? d : 2;
			h = h >= f? h : f;    // h = H(i,j) = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			h1 = h;              // save H(i,j) to h1 for the next column
			mj = m > h? mj : j;
			m = m > h? m : h;    // update the max score in this row
			h -= gapoe;
			e -= gape;
			d |= e > h? 1<<2 : 0;
			e = e > h? e : h;    // e = E(i+1,j)
			E[j] = e;            // save E(i+1,j) for the next row
			f -= gape;
			d |= f > h? 2<<4 : 0;
			f = f > h? f : h;    // f = F(i,j+1)
			z[i][j] = d;           // z[i,j] keeps h for the current cell and e/f for the next cell
			console.log(d);
		}
		H[query.length] = h1, E[query.length] = Number.NEGATIVE_INFINITY;
		if (m > max) max = m, end_i = i, end_j = mj;
	}
	score = H[query.length];

	// backtrack to recover the alignment/cigar
	function push_cigar(ci, op, len) {
		if (ci.length == 0 || op != (ci[ci.length-1]&0xf))
			ci.push(len<<4|op);
		else ci[ci.length-1] += len<<4;
	}
	var cigar = [], tmp, which = 0, i, k, start_i = 0;
	i = t.length - 1, k = (i + 1 < query.length? i + 1 : query.length) - 1; // (i,k) points to the last cell
	while (i >= 0 && k >= 0) {
		tmp = z[i][k];
		which = tmp >> (which << 1) & 3;
		if (which == 0 && tmp>>6) break;
		if (which == 0) which = tmp & 3;
		if (which == 0)      { push_cigar(cigar, 0, 1); --i, --k; } // match
		else if (which == 1) { push_cigar(cigar, 2, 1); --i; } // deletion
		else                 { push_cigar(cigar, 1, 1), --k; } // insertion
	}
	// add the first insertion or deletion
	if (i >= 0) push_cigar(cigar, 2, i + 1);
	if (k >= 0) push_cigar(cigar, 1, k + 1);
	cigar.reverse()
	return [score, start_i, cigar];
}

function printGaps(target, query, start, cigar)
{
	var oq = '', ot = '', mid = '', lq = 0, lt = start;
	for (var k = 0; k < cigar.length; ++k) {
		var op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 0) { // match
			oq += query.substr(lq, len);
			ot += target.substr(lt, len);
			lq += len, lt += len;
		} else if (op == 1) { // insertion
			oq += query.substr(lq, len);
			ot += Array(len+1).join("-");
			lq += len;
		} else if (op == 2) { // deletion
			oq += Array(len+1).join("-");
			ot += target.substr(lt, len);
			lt += len;
		} else if (op == 4) { // soft clip
			lq += len;
		}
	}
	var ut = ot.toUpperCase();
	var uq = oq.toUpperCase();
	for (var k = 0; k < ut.length; ++k)
		mid += ut.charAt(k) == uq.charAt(k)? '|' : ' ';
	return [ot, oq, mid];
}
