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
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFAlternativeAlleleCounter {

	private CloseableIterator<VariantContext>  iterator;
	ImmutableSet<IFilter> filters;
	private int counts = -1;

	public VCFAlternativeAlleleCounter(CloseableIterator<VariantContext> iterator, ImmutableSet<IFilter> filters) {
		this.iterator = iterator;
		this.filters = filters;
	}

	public int getCounts() {
		if (counts < 0)
			count();
		return counts;
	}

	private void count() {
		counts = 0;
		while (iterator.hasNext()) {
			VariantContext vc = iterator.next();
			for (IFilter iFilter : filters) {
				vc = iFilter.filter(vc);
			}
			if (vc != null)
				counts += vc.getAlternateAlleles().size();
		}
		// FIXME Really close it? test this!
		iterator.close();
	}
}
