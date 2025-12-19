use std::fmt::Display;

#[derive(Debug, Default)]
pub struct ChannelReporter {
    pub inp_num: usize,
    pub filter_by_input_filter: usize,
    pub filter_by_no_align: usize,
    pub filter_by_cs2ref_align_fail: usize,

    pub out_num: usize,
}

impl Display for ChannelReporter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "\
    inp_num                    :{},
    filter_by_input_filter     :{},
    filter_by_no_align         :{},
    filter_by_cs2ref_align_fail:{},
    oup_num                    :{:.4}% ({})
        ",
            self.inp_num,
            self.filter_by_input_filter,
            self.filter_by_no_align,
            self.filter_by_cs2ref_align_fail,
            (self.out_num as f32 / (self.inp_num as f32 + 1e-4)) * 100.,
            self.out_num,
        )
    }
}

#[derive(Debug, Default)]
pub struct SbrReporter {
    pub inp_num: usize,
    pub used_num: usize,
    pub filter_by_length: usize,
    pub filter_by_alignment: usize,
    pub fallback_num: usize,
    pub fallback_resuced_num: usize
}

impl Display for SbrReporter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "\
    inp_num:{},
    filter_by_length:{:.4}% ({})
    filter_by_alignment:{:.4}% ({})
    fallback_num:{:.4}% ({})
    fallback_rescued_num:{:.4}% ({})
        ",
            self.inp_num,
            (self.filter_by_length as f32 / (self.inp_num as f32 + 1e-4)) * 100.,
            self.filter_by_length,
            (self.filter_by_alignment as f32 / (self.inp_num as f32 + 1e-4)) * 100.,
            self.filter_by_alignment,
            (self.fallback_num as f32 / (self.inp_num as f32 + 1e-4)) * 100.,
            self.fallback_num,
            (self.fallback_resuced_num as f32 / (self.inp_num as f32 + 1e-4)) * 100.,
            self.fallback_resuced_num,

        )
    }
}

#[derive(Debug, Default)]
pub struct Reporter {
    pub channel_reporter: ChannelReporter,
    pub sbr_reporter: SbrReporter,
}

impl Display for Reporter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "\
ChannelYield:
    {}

SubreadsYield:
    {}
        ",
            self.channel_reporter, self.sbr_reporter
        )
    }
}


#[cfg(test)]
mod test {
    use super::Reporter;


    #[test]
    fn test_reporter() {

        let reporter = Reporter::default();
        println!("{}", reporter);

    }
}