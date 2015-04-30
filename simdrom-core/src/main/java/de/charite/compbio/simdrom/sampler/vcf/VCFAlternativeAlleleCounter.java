/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class VCFAlternativeAlleleCounter {

	private VCFFileReader parser;
	private int counts = -1;

	public VCFAlternativeAlleleCounter(String filePath) {
		this.parser = new VCFFileReader(new File(filePath), false);
	}

	public int getCounts() {
		if (counts < 0)
			count();
		return counts;
	}

	private void count() {
		counts = 0;
		CloseableIterator<VariantContext> iterator = parser.iterator();
		while (iterator.hasNext()) {
			counts += iterator.next().getAlternateAlleles().size();
			iterator.next();
		}
		parser.close();
	}
}
