package de.charite.compbio.simdrom;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

import java.io.IOException;

import org.apache.commons.cli.ParseException;

import de.charite.compbio.simdrom.cli.SIMdromSetting;
import de.charite.compbio.simdrom.io.writer.VCFTSVWriter;
import de.charite.compbio.simdrom.sampler.SpikeIn;
import de.charite.compbio.simdrom.sampler.vcf.VCFAlternativeAlleleCounter;
import de.charite.compbio.simdrom.sampler.vcf.VCFRandomSampleSelecter;
import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class Main {

	public static void main(String[] args) throws ParseException, IOException {

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

		// writer
		VariantContextWriter writer = new VariantContextWriterBuilder().setOutputVCFStream(System.out)
				.unsetOption(Options.INDEX_ON_THE_FLY).build();

		// log
		boolean log = SIMdromSetting.SPLIKE_IN_LOGFILE != null;

		SpikeIn spikein = new SpikeIn(backgroundSampler, mutationSampler, log);

		// header
		writer.writeHeader(spikein.getVCFHeader());

		// spike in and write out
		while (spikein.hasNext()) {
			writer.add(spikein.next());
		}

		// write log
		if (log) {
			VCFTSVWriter logWriter = new VCFTSVWriter(SIMdromSetting.SPLIKE_IN_LOGFILE);
			boolean header = false;
			for (VariantContext vc : spikein.getVcLogs()) {
				if (!header) {
					logWriter.writeHeader(vc);
					header = true;
				}
				logWriter.add(vc);
			}
			logWriter.close();
		}

		// close properly
		writer.close();
		spikein.close();
		// and exit properly
		System.exit(0);
	}
}
