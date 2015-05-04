package de.charite.compbio.simdrom.io.writer;

import htsjdk.variant.variantcontext.VariantContext;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import com.google.common.collect.ImmutableList;

public class VCFTSVWriter implements Closeable {

	private final CSVFormat format = CSVFormat.newFormat('\t').withQuote(null).withRecordSeparator("\r\n")
			.withIgnoreSurroundingSpaces(true);
	private CSVPrinter printer;
	private ImmutableList<String> header;
	private static final ImmutableList<String> header_start = ImmutableList.of("#CHROM", "POS", "ID", "REF", "ALT",
			"QUAL", "FILTER");

	public VCFTSVWriter(String file) throws IOException {
		this.printer = new CSVPrinter(new BufferedWriter(new FileWriter(file)), format);
	}

	public void setHeader(VariantContext vc) {
		this.header = ImmutableList.<String> builder().addAll(header_start)
				.addAll(vc.getCommonInfo().getAttributes().keySet()).build();
	}

	public void writeHeader(VariantContext vc) throws IOException {
		if (header == null)
			setHeader(vc);
		printer.printRecord(this.header);
	}

	public void add(VariantContext vc) throws IOException {
		List<Object> record = new ArrayList<Object>();
		record.add(vc.getContig());
		record.add(vc.getStart());
		record.add(vc.getID());
		record.add(vc.getReference());
		record.add(vc.getAlternateAlleles());
		record.add(vc.getPhredScaledQual());
		record.add(vc.getFilters());
		for (String header : getHeader()) {
			if (!header_start.contains(header)) {
				if (vc.getCommonInfo().hasAttribute(header)) {
					record.add(vc.getCommonInfo().getAttribute(header));
				}
			}
		}
		this.printer.printRecord(record);
	}

	public ImmutableList<String> getHeader() {
		if (header == null)
			this.header = ImmutableList.<String> builder().addAll(header_start).build();
		return header;
	}

	@Override
	public void close() throws IOException {
		printer.close();

	}

}
