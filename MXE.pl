use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Transcript;
use Bio::Seq;
use Bio::SeqIO;

use Data::Dumper;

print "Starting...\n";

my $registry = 'Bio::EnsEMBL::Registry';


$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous',
	-species => 'mus musculus'
);

my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );

print "Registry loaded.\n";

my $totalMXE = 0; #number of total MXE pairs
my $noValidTranslations = 0; #number of MXE pairs which result in 0 translations without stop codons
my $atLeastOneValidTranslation = 0; #number of MXE pairs which have at least one translation without stop codons
my $atLeastOneValidTranslationPair = 0; #number of MXE pairs which have at least one valid pair of translations in the same frame without stop codons
my $exactlyOneValidTranslationPair = 0; #number of MXE pairs which have exactly one valid pair of translations in the same frame, and no other translations

my $strandIsCorrect = 0;
my $strandIncorrect = 0;


my $out = Bio::SeqIO->new(-file => ">fasta.txt" , '-format' => 'fasta');


while(<>)
{

	my $linestring = $_;

	my @values = split("\t", $linestring);

	my $column1 = $values[0];
	
	next if ($column1 eq "ID");

	my $geneID = $values[1];
	my $geneSymbol = $values[2];
	$geneID =~ tr/"//d;
	$geneSymbol =~ tr/"//d;
	my $chromosome = $values[3];
	my $strand = $values[4];
	
	my $stStart = $values[5];
	my $stEnd = $values[6];
	my $ndStart = $values[7];
	my $ndEnd = $values[8];
	my $upstreamStart  = $values[9];
	my $upstreamEnd  = $values[10];
	my $downstreamStart  = $values[11];
	my $downstreamEnd  = $values[12];

	my $slice1 = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $stStart + 1, $stEnd );
	my $slice2 = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $ndStart + 1, $ndEnd  );
	my $upstreamSlice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $upstreamStart + 1, $upstreamEnd  );
	my $downstreamSlice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $downstreamStart + 1, $downstreamEnd );

	#### STEP 1: Check to see if UTRs overlap either the upstream or the downstream slices

	
	my @transcripts = @{$upstreamSlice->get_all_Transcripts()};
	
	my $upsteamContainsFivePrime = 0;
	my $downstreamContainsThreePrime = 0;
	
	my $fivePrimeEnd = -1;
	my $threePrimeStart = 9007199254740992; #max int
	
	foreach $t (@transcripts){
	
		bless $t, "Bio::EnsEMBL::Transcript";
		
		my $five_prime  = $t->five_prime_utr_Feature
                 or next;
		
		$five_prime = $five_prime->transform('chromosome');
		
		if($five_prime->end > $upstreamStart and $five_prime->end < $upstreamEnd)
		{
			$upsteamContainsFivePrime = 1;
			$fivePrimeEnd = $fivePrimeEnd >= $five_prime->end ? $fivePrimeEnd : $five_prime->end;
		}
	}
	
	foreach $t (@transcripts){
		
		my $three_prime  = $t->three_prime_utr_Feature
                 or next;
		
		$three_prime = $three_prime->transform('chromosome');
		
		if($three_prime->start < $downstreamEnd and $three_prime->start > $downstreamStart)
		{
			$downstreamContainsThreePrime = 1;
			$threePrimeStart = $threePrimeStart <= $three_prime->start ? $threePrimeStart : $three_prime->start;
		}
	}
	
	#### STEP 2: If UTRs overlap either the upstream or the downstream slices, trim these slices to remove the overlap
	
	if($fivePrimeEnd > 1)
	{
		$upstreamSlice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $fivePrimeEnd + 1, $upstreamEnd  );
		print "\nupstream slice reduced\n";
	}
	
	if($downstreamContainsThreePrime)
	{
		$downstreamSlice = $slice_adaptor->fetch_by_region( 'chromosome', $chromosome, $downstreamStart + 1, $threePrimeStart - 1 );
		print "\ndownstream slice reduced\n";
	}
	
	
	
	#### STEP 3: Generate sequences from the (trimmed) slices
	
	my $sequence1 = $upstreamSlice->seq().$slice1->seq().$downstreamSlice->seq();
	my $sequence2 = $upstreamSlice->seq().$slice2->seq().$downstreamSlice->seq();

	#### STEP 4: Validate the sequences
	
	print "\n\nValid Protein Sequences:";

	my $seq1_obj = Bio::Seq->new(-seq => $sequence1, -alphabet => 'dna' );
	my $seq2_obj = Bio::Seq->new(-seq => $sequence2, -alphabet => 'dna' );
	$seq1_obj->id($geneSymbol."-"."MXE"."-1"."-".$column1."-".$geneID);
	$seq2_obj->id($geneSymbol."-"."MXE"."-2"."-".$column1."-".$geneID);
	$seq1_obj->accession_number($geneID);
	$seq2_obj->accession_number($geneID);
	
	my @seqs = Bio::SeqUtils->translate_6frames($seq1_obj);
	my @seqs2 = Bio::SeqUtils->translate_6frames($seq2_obj);
	
	my $validTranslations = 0;
	my $validTranslationPairs = 0;
	for(my $i = 0; $i < @seqs; $i++)
	{
		my $seq1 = $seqs[$i];
		my $seq2 = $seqs2[$i];
		
		next if ($strand eq "-" and $seq1->id =~ m/[0-2]F/);
		next if ($strand eq "+" and $seq1->id =~ m/[0-2]R/);
		
		my $seq1good = not $seq1->seq =~ m/\*/;
		my $seq2good = not $seq2->seq =~ m/\*/;
		
		if($seq1good)
		{
			$validTranslations++;
			print "\n\n", $seq1->id, " ", $seq1->seq;
		}
		if($seq2good)
		{
			$validTranslations++;
			print "\n\n", $seq2->id, " ", $seq2->seq;
		}
		if($seq1good and $seq2good)
		{
		
			# if ($strand eq "-" and $seq1->id =~ m/[0-2]R/)
			# {
				# $strandIsCorrect++;
			# }
			# if ($strand eq "+" and $seq1->id =~ m/[0-2]F/)
			# {
				# $strandIsCorrect++;
			# }
			
			# if ($strand eq "-" and $seq1->id =~ m/[0-2]F/)
			# {
				# $strandIncorrect++;
			# }
			# if ($strand eq "+" and $seq1->id =~ m/[0-2]R/)
			# {
				# $strandIncorrect++;
			# }
		
			$out->write_seq($seq1);
			$out->write_seq($seq2);
			$validTranslationPairs++;
		}
	}
	
	$totalMXE++;
	if($validTranslations)
	{
		$atLeastOneValidTranslation++;
	}
	else
	{
		print " None";
		$noValidTranslations++;
		
	}

	if($validTranslationPairs)
	{
		$atLeastOneValidTranslationPair++;
	}

	if($validTranslationPairs == 1 and $validTranslations == 2)
	{
		$exactlyOneValidTranslationPair++;		
	}
	
	if($totalMXE % 10 == 0)
	{
		print "\n\n$totalMXE processed so far"
	}
	
}

print "\n\ntotalMXE $totalMXE\n";
print "noValidTranslations $noValidTranslations\n";
print "atLeastOneValidTranslation $atLeastOneValidTranslation\n";
print "atLeastOneValidTranslationPair $atLeastOneValidTranslationPair\n";
print "exactlyOneValidTranslationPair $exactlyOneValidTranslationPair\n";

#print "strand correct $strandIsCorrect\n";
#print "strand incorrect $strandIncorrect";