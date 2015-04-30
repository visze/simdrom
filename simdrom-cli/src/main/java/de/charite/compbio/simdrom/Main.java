package de.charite.compbio.simdrom;

import java.util.LinkedHashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.apache.commons.cli.ParseException;

import de.charite.compbio.simdrom.cli.SIMdromSetting;
import de.charite.compbio.simdrom.sampler.vcf.VCFAlternativeAlleleCounter;
import de.charite.compbio.simdrom.sampler.vcf.VCFRandomSampleSelecter;
import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class Main {

	public static void main(String[] args) throws ParseException {

		SIMdromSetting.parse(args);

		VCFSampler backgroundSampler = new VCFSampler(SIMdromSetting.BACKGROUND_VCF);

		backgroundSampler.setProbability(SIMdromSetting.BACKGROUND_PROBABILITY);

		if (SIMdromSetting.BACKGROUND_PROBABILITY > 1.0) {
			VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(SIMdromSetting.BACKGROUND_VCF);
			backgroundSampler.setCounts(counter.getCounts());
		}

		if (SIMdromSetting.ONLY_ONE_SAMPLE) {
			VCFRandomSampleSelecter selecter = new VCFRandomSampleSelecter(SIMdromSetting.BACKGROUND_VCF);
			backgroundSampler.setSample(selecter.getSample());
		}

		if (SIMdromSetting.ALLELE_FREQUENCY_IDENTIFIER != null) {
			backgroundSampler.setAFIdentifier(SIMdromSetting.ALLELE_FREQUENCY_IDENTIFIER);
		}

		VCFSampler mutationSampler = null;
		if (SIMdromSetting.MUTATIONS_VCF != null) {
			mutationSampler = new VCFSampler(SIMdromSetting.MUTATIONS_VCF);

			mutationSampler.setProbability(SIMdromSetting.MUTATIONS_PROBABILITY);

			if (SIMdromSetting.MUTATIONS_PROBABILITY > 1.0) {
				VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(SIMdromSetting.MUTATIONS_VCF);
				mutationSampler.setCounts(counter.getCounts());
			}
		}

		VariantContextWriter writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out)
				.unsetOption(Options.INDEX_ON_THE_FLY).build();
		
//		header
		if (mutationSampler != null) {
			Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>();
			metaData.addAll(backgroundSampler.getFileHeader().getMetaDataInInputOrder());
			metaData.addAll(mutationSampler.getFileHeader().getMetaDataInInputOrder());
			writer.writeHeader(new VCFHeader(metaData));
			
		} else 
			writer.writeHeader(backgroundSampler.getFileHeader());
		
		VariantContext backgroundVC = null;
		VariantContext mutationsVC = null;
		
		if (backgroundSampler.hasNext())
			backgroundVC = backgroundSampler.next();
		if ((mutationSampler != null && mutationSampler.hasNext()))
			mutationsVC = mutationSampler.next();
		while (backgroundVC != null || mutationsVC != null) {
			boolean backgroundSelection = true;
			
			if (mutationsVC == null)
				writer.add(backgroundVC);
			else if (backgroundVC == null) {
				writer.add(mutationsVC);
				backgroundSelection = false;
			} else if (backgroundVC.getContig().equals(mutationsVC.getContig())) {
				if (backgroundVC.getStart() <= mutationsVC.getStart())
					writer.add(backgroundVC);
				else {
					writer.add(mutationsVC);
					backgroundSelection = false;
				}
			} else {
				writer.add(backgroundVC);
			}
			
			if (backgroundSampler.hasNext() && backgroundSelection)
				backgroundVC = backgroundSampler.next();
			if ((mutationSampler != null && mutationSampler.hasNext()) && !backgroundSelection)
				mutationsVC = mutationSampler.next();
		}
		writer.close();
		backgroundSampler.close();
		if (mutationSampler != null)
			mutationSampler.close();
	}
}
