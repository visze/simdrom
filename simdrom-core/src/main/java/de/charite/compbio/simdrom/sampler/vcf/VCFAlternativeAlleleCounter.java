/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.filter.IFilter;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class VCFAlternativeAlleleCounter {

	private VCFFileReader parser;
	 ImmutableSet<IFilter> filters;
	private int counts = -1;

	public VCFAlternativeAlleleCounter(String filePath, ImmutableSet<IFilter> filters) {
		this.parser = new VCFFileReader(new File(filePath), false);
		this.filters = filters;
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
			VariantContext vc = iterator.next();
			for (IFilter iFilter : filters) {
				vc = iFilter.filter(vc);
			}
			if (vc != null)
				counts += vc.getAlternateAlleles().size();
		}
		parser.close();
	}
}
