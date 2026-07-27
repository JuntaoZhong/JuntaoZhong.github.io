window.onload = initialize;

function initialize() {
	// Set up the "help" button's event handler.
	document.getElementById("information_sign_button").onclick = informationSign;
}

function informationSign(){
	var b = document.getElementById("information_sign_button");
	var c = document.getElementById("information_sign_content"); 
	if (b.value == "What is DNA sequence alignment???") {
		window.scrollTo(0,0);
		c.style.display = "block";
		b.value = "Hide Help";
	}else{
		c.style.display = "none";
		b.value = "What is DNA sequence alignment???";
	}
}

function do_align(){
	var time_start = new Date().getTime();

	var target = document.getElementById('target').value.replace(/\s+/g, '').replace('\n', '');
	var query  = document.getElementById('query').value.replace(/\s+/g, '').replace('\n', '');
	var matchScore   = parseInt(document.getElementById('match').value);
	var mismatchScore  = parseInt(document.getElementById('mismatch').value);
	var gapOpen = parseInt(document.getElementById('gapOpen').value);
	var gapExtent = parseInt(document.getElementById('gapExtent').value);

	var rst = align(target, query, [matchScore, mismatchScore], [gapOpen, gapExtent]);
	var str = 'alignment score: ' + rst[0] + '\n';
	str += 'alignment:\n';
	var fmt = printGaps(target, query, rst[1], rst[2]);

	var linelen = 100, n_lines = 8;
	for (var l = 0; l < fmt[0].length; l += linelen) {
		str += fmt[0].substr(l, linelen) + '\n';
		str += fmt[2].substr(l, linelen) + '\n';
		str += fmt[1].substr(l, linelen) + '\n\n';
	}

	document.getElementById('out').value = str;
	document.getElementById('out').rows = n_lines;

	var elapse = (new Date().getTime() - time_start) / 1000.0;
	document.getElementById('runtime').innerHTML = "in " + elapse.toFixed(3) + "s";
}

var humanDNA=`ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCG\
CACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACC\
AAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGC\
AAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCG\
CTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAG\
CTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCT\
GCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAA\
TACCGTTAA`

var fishDNA=`ATGAGTCTCACTGCCAAGGACAAGGAAACAGTCAAAGCCTTCTGGGCTAAAGTGGCTCCC\
AAGGCTGAAGACATTGGCCAGGATGCTCTGTCCAGGATGCTGGCGGTTTACCCACAGACC\
AAGACCTACTTCTCCCACTGGAAGGACATGAGTGCCGGCTCTGCTCCAGTGAAGAAGCAC\
GGAGCTACGGTGATGGGTGGAGTAGCTGATGCTGTGACCAAAATCGATGATCTGACCTCA\
GGTCTCCTGAGCCTGAGTGAGCTGCATGCTTTCACTCTTAGAGTGGACCCTGCCAACTTC\
AAGATCCTGGCACACAACATCCTTGTGGTCTTCGCCATCAAGTTTCCCACCGACTTCACC\
CCTGAGGTCCATGTGTCTGTGGACAAGTTCTTGGCTGCTCTGGCCCGAGCCCTCTCCGAG\
AAGTACAGATAA`

