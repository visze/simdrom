package de.charite.compbio.simdrom.interval;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

public class SAMFileHeaderBuilder {

	public static SAMFileHeader build() {

		SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("1", 249250621));
		header.addSequence(new SAMSequenceRecord("2", 243199373));
		header.addSequence(new SAMSequenceRecord("3", 198022430));
		header.addSequence(new SAMSequenceRecord("4", 191154276));
		header.addSequence(new SAMSequenceRecord("5", 180915260));
		header.addSequence(new SAMSequenceRecord("6", 171115067));
		header.addSequence(new SAMSequenceRecord("7", 159138663));
		header.addSequence(new SAMSequenceRecord("8", 146364022));
		header.addSequence(new SAMSequenceRecord("9", 141213431));
		header.addSequence(new SAMSequenceRecord("10", 135534747));
		header.addSequence(new SAMSequenceRecord("11", 135006516));
		header.addSequence(new SAMSequenceRecord("12", 133851895));
		header.addSequence(new SAMSequenceRecord("13", 115169878));
		header.addSequence(new SAMSequenceRecord("14", 107349540));
		header.addSequence(new SAMSequenceRecord("15", 102531392));
		header.addSequence(new SAMSequenceRecord("16", 90354753));
		header.addSequence(new SAMSequenceRecord("17", 81195210));
		header.addSequence(new SAMSequenceRecord("18", 78077248));
		header.addSequence(new SAMSequenceRecord("19", 59128983));
		header.addSequence(new SAMSequenceRecord("20", 63025520));
		header.addSequence(new SAMSequenceRecord("21", 48129895));
		header.addSequence(new SAMSequenceRecord("22", 51304566));
		header.addSequence(new SAMSequenceRecord("MT", 16569));
		header.addSequence(new SAMSequenceRecord("M", 16569));
		header.addSequence(new SAMSequenceRecord("X", 155270560));
		header.addSequence(new SAMSequenceRecord("Y", 59373566));
		header.addSequence(new SAMSequenceRecord("chr1", 249250621));
		header.addSequence(new SAMSequenceRecord("chr2", 243199373));
		header.addSequence(new SAMSequenceRecord("chr3", 198022430));
		header.addSequence(new SAMSequenceRecord("chr4", 191154276));
		header.addSequence(new SAMSequenceRecord("chr5", 180915260));
		header.addSequence(new SAMSequenceRecord("chr6", 171115067));
		header.addSequence(new SAMSequenceRecord("chr7", 159138663));
		header.addSequence(new SAMSequenceRecord("chr8", 146364022));
		header.addSequence(new SAMSequenceRecord("chr9", 141213431));
		header.addSequence(new SAMSequenceRecord("chr10", 135534747));
		header.addSequence(new SAMSequenceRecord("chr11", 135006516));
		header.addSequence(new SAMSequenceRecord("chr12", 133851895));
		header.addSequence(new SAMSequenceRecord("chr13", 115169878));
		header.addSequence(new SAMSequenceRecord("chr14", 107349540));
		header.addSequence(new SAMSequenceRecord("chr15", 102531392));
		header.addSequence(new SAMSequenceRecord("chr16", 90354753));
		header.addSequence(new SAMSequenceRecord("chr17", 81195210));
		header.addSequence(new SAMSequenceRecord("chr18", 78077248));
		header.addSequence(new SAMSequenceRecord("chr19", 59128983));
		header.addSequence(new SAMSequenceRecord("chr20", 63025520));
		header.addSequence(new SAMSequenceRecord("chr21", 48129895));
		header.addSequence(new SAMSequenceRecord("chr22", 51304566));
		header.addSequence(new SAMSequenceRecord("chrMT", 16569));
		header.addSequence(new SAMSequenceRecord("chrM", 16569));
		header.addSequence(new SAMSequenceRecord("chrX", 155270560));
		header.addSequence(new SAMSequenceRecord("chrY", 59373566));

		return header;
	}
}
