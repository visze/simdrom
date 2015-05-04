package de.charite.compbio.simdrom.io.writer;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFEncoder;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
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
		List<String> record = new ArrayList<String>();
		record.add(vc.getContig());
		record.add(String.valueOf(vc.getStart()));
		record.add(vc.getID());
		// REF
		record.add(vc.getReference().getDisplayString());
		// ALT
		final StringBuilder stringBuilder = new StringBuilder();
		if (vc.isVariant()) {
			Allele altAllele = vc.getAlternateAllele(0);
			String alt = altAllele.getDisplayString();
			stringBuilder.append(alt);

			for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
				altAllele = vc.getAlternateAllele(i);
				alt = altAllele.getDisplayString();
				stringBuilder.append(",");
				stringBuilder.append(alt);
			}
			record.add(stringBuilder.toString());
		} else {
			record.add(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
		}
		// QUAL
		if (!vc.hasLog10PError())
			record.add(VCFConstants.MISSING_VALUE_v4);
		else
			stringBuilder.append(formatQualValue(vc.getPhredScaledQual()));
		// FILTER
		record.add(formatVCFField(vc.getFilters()));
		// INFO
		for (String header : getHeader()) {
			if (!header_start.contains(header)) {
				if (vc.getCommonInfo().hasAttribute(header)) {
					record.add(formatVCFField(vc.getCommonInfo().getAttribute(header)));
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

	private String formatVCFField(final Object val) {
		final String result;
		if (val == null)
			result = VCFConstants.MISSING_VALUE_v4;
		else if (val instanceof Double)
			result = VCFEncoder.formatVCFDouble((Double) val);
		else if (val instanceof Boolean)
			result = (Boolean) val ? "" : null; // empty string for true, null
												// for false
		else if (val instanceof List) {
			result = formatVCFField(((List) val).toArray());
		} else if (val.getClass().isArray()) {
			final int length = Array.getLength(val);
			if (length == 0)
				return formatVCFField(null);
			final StringBuilder sb = new StringBuilder(formatVCFField(Array.get(val, 0)));
			for (int i = 1; i < length; i++) {
				sb.append(",");
				sb.append(formatVCFField(Array.get(val, i)));
			}
			result = sb.toString();
		} else
			result = val.toString();

		return result;
	}

	private static final String QUAL_FORMAT_STRING = "%.2f";
	private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

	private String formatQualValue(final double qual) {
		String s = String.format(QUAL_FORMAT_STRING, qual);
		if (s.endsWith(QUAL_FORMAT_EXTENSION_TO_TRIM))
			s = s.substring(0, s.length() - QUAL_FORMAT_EXTENSION_TO_TRIM.length());
		return s;
	}

}
